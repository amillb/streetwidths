#!/usr/bin/python

"""
Process a county streets and parcels data
to calculate land values
"""

import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np
from osgeo import ogr
import os
import seaborn as sns
import shutil
import time
from pgMapMatch import tools

figsizes = {'full': (6.5,5), 'short':(6.5,3),'halfwidth':(3.25,3)}
sqm_per_acre, m2ft = 4046.86, 3.28084
c5s = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00','k']
c2s = ['#7570b3','#1b9e77']

# projections are State Plane (NAD83(NSRS2007)) meters, code is EPSG
projections = {'alameda_ca':3493, 'bexar_tx':3673, 'cook_il':3528, 'dallas_tx':3669, 'harris_tx':3673, 
               'hennepin_mn':3596, 'king_wa':3689, 'kings_ny':3627, 'los_angeles_ca':3497, 'maricopa_az':3478,
               'miami_dade_fl':3511, 'middlesex_ma':3585, 'orange_ca':3499, 'queens_ny':3627, 'riverside_ca':3499,
            'san_bernardino_ca':3497, 'san_diego_ca':3499, 'santa_clara_ca':3493, 'shelby_tn':3661, 'tarrant_tx':3669}
fipslookup = {'alameda_ca':'06001', 'bexar_tx':'48029', 'cook_il':'17031', 'dallas_tx':'48113', 'harris_tx':'48201', 
              'hennepin_mn':'27053', 'king_wa':'53033', 'kings_ny':'36047', 'los_angeles_ca':'06037', 'maricopa_az':'04013',
            'miami_dade_fl':'12086', 'middlesex_ma':'25017', 'orange_ca':'06059', 'queens_ny':'36081', 'riverside_ca':'06065',
            'san_bernardino_ca':'06071', 'san_diego_ca':'06073', 'santa_clara_ca':'06085', 'shelby_tn':'47157','tarrant_tx':'48439'}
countylookup = {fips:county for county, fips in fipslookup.items()}
def format_county(cname):
    cname = cname.replace('_',' ').title()
    return cname[:-3]+','+cname[-3:].upper() # format state name
countylookup_formatted = {fips:format_county(county) for county, fips in fipslookup.items()}

def definePaths():
    basepath = '/Users/adammb/Documents/GDrive/Research/Street widths/'
    datapath = '/Volumes/Data/Streets/'
    paths = {'data':datapath+'Raw_data/',
             'atlas':basepath+'Data/',
             'code':basepath+'Code/',
             'scratch':datapath+'Scratch/',
             'working':basepath+'Working/',
             'graphics':basepath+'Figures/',
             'website':basepath+'Streets_website/',}
    return paths
paths = definePaths()

class logger():
    """
    Keep a log of progress/status messages.  Could be rewritten to use python's logging module?
    Optionally displays time between messages.
    Optionally prints to stdout as well.

    Passing topic=None skips any file writing, but still (if verbose = True) writes to stdout.
    """
    def __init__(self,topic, verbose=True, timing=False, prefix=None, overwrite=False):
        """
        If topic  is None, then no file is written.

        overwrite = True will delete/recreate the log file instead of appending to it
        """
        self.prefix='' if prefix is None else prefix+': '
        assert timing in [True,False]
        self.timing=timing
        logdir=paths['working']+'logs/'
        if not os.path.exists(logdir):
            os.mkdir(logdir)
        if topic is None:
            self.fout=open(os.devnull,"w")
            topic='null logger: no file recording'
        else:
            fname=logdir+topic+'.log'
            self.fout=open(fname, ['a','w'][int(overwrite)])
        self.verbose=verbose
        self.startime=time.time()
        self.latesttime=self.startime
        self.write('---------- NEW LOGGER: '+self.prefix+topic+':   '+time.strftime("%Y:%m:%d:%H:%M:%S")+' ------------\n')
    def prependTime(self,sss):
        if not self.timing: return sss
        nowt=time.time()
        dtime=nowt-self.latesttime
        etime=nowt-self.startime
        self.latesttime=nowt
        return  '['+time.strftime("%Y:%m:%d:%H:%M:%S", time.localtime(nowt))+'][At '+('%5d s'%etime if etime<61 else '%5.1f m'%(etime*1.0/60) )+'|'+ 'dt='+('%7.2g s'%dtime if dtime<61 else ' %7.2f m'%(dtime*1.0/60) )   + '] '+sss

    def write(self,sss):
        outs=self.prependTime(self.prefix+sss)
        if self.verbose:
            print(outs)
        self.fout.write(outs+'\n')
    def close(self):
        self.fout.close()

class dataLoader():
    """Uploads parcel files and streets into Postgres"""

    def __init__(self, forceUpdate=False):
        self.db = tools.dbConnection(user='adammb', db='streetwidths', host='localhost', schema='rawdata', requirePassword=False)
        self.forceUpdate = forceUpdate
   
        # some counties have an extra csv or dbf data file with year built
        # format is {fips: [filename, joinfield_shp, joinfield_extra yrbuiltfield, separator]}
        self.extradata = {'06037': ['Assessor_Parcels_Data_-_2017.csv',['ain'],['AIN'],'YearBuilt',','],
                          '06065': ['Assessor_PropertyCharacteristics_Table.csv',['APN'],['ASMT'],'YEAR_CONSTRUCTED',','],
                          '48113': ['res_detail.csv', ['Acct'],['ACCOUNT_NUM'],'YR_BUILT',','],
                          '48439': ['PropertyData.txt', ['taxpin'],['GIS_Link'],'Year_Built','|'],
                          '53033': ['EXTR_ResBldg.csv', ['MAJOR','MINOR'], ['Major','Minor'],'YrBuilt',',']}
        # parcel identifiers for right of way. first is original field name (may be capitalized). second is query
        self.rowParcel = {'48201':['parcel_typ','parcel_typ=6'], '06065':['APN',"apn='RW'"],'48439':['PARCELTYPE','parceltype=3'] }

    def loadParcelsandRoads(self):
        """Loads parcel file to PostGres"""
        tables = self.db.list_tables()
        for county in projections:
            fips = fipslookup[county]
            proj = projections[county]

            if 'parcels_'+fips in tables and not self.forceUpdate:
                continue

            shp = [ff for ff in os.listdir(paths['data']+'Parcel_data/'+county) if ff.endswith('.shp')]
            if len(shp)>1: 
                raise Exception('Multiple shp files found for {}'.format(county))
            shpFn = paths['data'] + 'Parcel_data/{}/{}'.format(county, shp[0])
            
            # get field names. Need at least one field (e.g. FID), plus year built and join field if available
            ds = ogr.Open(shpFn)
            lyr = ds.GetLayer()
            field_names = [field.name for field in lyr.schema if field.name.lower() in ['apn','fid','yearbuilt','yr_built','year_built','yrblt','yr_blt','year_constructed']]
            if field_names == []: # need at least one
                field_names = [field.name for field in lyr.schema][0:1]
            if fips in self.extradata:
                field_names += [fn for fn in self.extradata[fips][1] if fn not in field_names]
            if fips in self.rowParcel and self.rowParcel[fips][0] not in field_names:
                field_names += [self.rowParcel[fips][0]] 

            field_names = ' '.join(field_names)
            #print(' '.join([field.name for field in lyr.schema]))
            ds = None # to close it

            print('Loading {}: {}'.format(county, fips))
            cmd = '''ogr2ogr -f "PostgreSQL" "PG:dbname=streetwidths user=adammb" {} -overwrite -nln rawdata.parcels_{} -select "{}" -nlt PROMOTE_TO_MULTI -lco GEOMETRY_NAME=geom -t_srs EPSG:{}'''.format(shpFn, fips, field_names, proj)
            assert os.system(cmd) == 0
             
            if fips in self.rowParcel:  # delete ROW parcels
                self.db.execute('DELETE FROM rawdata.parcels_{} WHERE {};'.format(fips, self.rowParcel[fips][1]))

            self.db.execute('CREATE INDEX pcl_{}_spat_idx ON rawdata.parcels_{} USING gist (geom);'.format(fips, fips))

            # adjust geometry of invalid parcels
            self.db.execute('''DELETE FROM rawdata.parcels_{} WHERE NOT ST_IsValid(geom) AND ST_Area(geom)<1;'''.format(fips))
            print('Dropped {} invalid parcels with zero area'.format(self.db.cursor.rowcount))
            # making it valid creates a line string and a polygon
            self.db.execute('''UPDATE rawdata.parcels_{} SET geom=ST_CollectionExtract(ST_MakeValid(geom), 3) -- get only polygons, ditch the line
                             WHERE NOT ST_IsValid(geom) 
                             AND ST_IsValidReason(geom) LIKE 'IllegalArgumentException: Invalid number of points in LinearRing%';'''.format(fips))
            print('Updated {} parcels with invalid number of points'.format(self.db.cursor.rowcount))
            # Buffer self intersections. See https://postgis.net/workshops/postgis-intro/validity.html
            self.db.execute('''UPDATE rawdata.parcels_{} SET geom=ST_Multi(ST_Buffer(geom,0)) 
                             WHERE NOT ST_IsValid(geom) 
                             AND ST_IsValidReason(geom) LIKE '%Self-intersection%';'''.format(fips))
            print('Updated {} invalid parcels with self-intersection'.format(self.db.cursor.rowcount))
            #self.db.execute('''UPDATE rawdata.parcels_{} SET geom=ST_Multi(ST_ConvexHull(geom)) WHERE NOT ST_IsValid(geom);'''.format(fips))
            #print('Updated {} invalid parcels'.format(self.db.cursor.rowcount))
            assert self.db.execfetch('''SELECT COUNT(*) FROM rawdata.parcels_{} WHERE NOT ST_IsValid(geom);'''.format(fips))[0][0]==0

            if fips == '25017': # middlesex is whole state. Restrict to that county
                cmd = '''CREATE TABLE rawdata.parcels_25017_tmp AS
                         SELECT t1.*
                         FROM rawdata.parcels_25017 t1, rawdata.tracts t2
                         WHERE ST_Intersects(t1.geom, ST_Transform(t2.geom,{proj})) AND t2.statefp='25' AND t2.countyfp='017';
                         DROP TABLE rawdata.parcels_25017;
                         ALTER TABLE rawdata.parcels_25017_tmp RENAME TO parcels_25017;
                         CREATE INDEX pcl_25017_spat_idx ON rawdata.parcels_25017 USING gist (geom);'''.format(proj=proj)
                self.db.execute(cmd)

            # add year built from separate file, if exists
            if fips in self.extradata:
                print('Loading extra data for {}'.format(county))
                joinFlds_shp, joinFlds_extra, yrfield, sep = self.extradata[fips][1:]
                joinFlds_shp = [jf.lower() for jf in joinFlds_shp] # will have been forced to lower case
                df = pd.read_csv(paths['data']+'Parcel_data/'+county+'/'+self.extradata[fips][0], sep=sep, usecols=joinFlds_extra+[yrfield], dtype={jf:str for jf in joinFlds_extra})
                df = df.groupby(joinFlds_extra)[yrfield].min().reset_index()  # take minimum if multiple for same parcel
                if joinFlds_shp!=joinFlds_extra:
                    df.rename(columns={jf1:jf2 for jf1, jf2 in zip(joinFlds_extra, joinFlds_shp)}, inplace=True)
                df.columns = [cc.lower() for cc in df.columns]
       
                self.db.df2db(df,'tmpjoin',index=False)
                self.db.execute('CREATE INDEX tmpjoin_idx ON rawdata.tmpjoin ({});'.format(','.join(joinFlds_shp)))
                self.db.execute('CREATE INDEX tmpjoin2_idx ON rawdata.parcels_{} ({});'.format(fips, ','.join(joinFlds_shp)))
                self.db.merge_table_into_table('tmpjoin', 'parcels_{}'.format(fips), joinFlds_shp)
                self.db.execute('DROP TABLE rawdata.tmpjoin;')
                self.db.execute('DROP INDEX tmpjoin2_idx;')                

            # rename to standard column name for yearbuilt
            cols = self.db.list_columns_in_table('parcels_{}'.format(fips))
            for col in ['year_built','yr_built','yrblt','yr_blt','year_constructed']:
                if col in cols:
                    self.db.execute('ALTER TABLE rawdata.parcels_{} RENAME {} to yearbuilt;'.format(fips, col))

            if 0:   # loads the census roads, but we don't use them 
                cmd = '''ogr2ogr -f "PostgreSQL" "PG:dbname=streetwidths user=adammb" {}Roads/tl_2019_{}_roads.shp -overwrite -nln rawdata.roads_{} -nlt PROMOTE_TO_MULTI -lco GEOMETRY_NAME=geom -t_srs EPSG:{}'''.format(paths['data'], fips, fips, proj)
                assert os.system(cmd) == 0
                self.db.execute('CREATE INDEX roads_{}_spat_idx ON rawdata.roads_{} USING gist (geom);'.format(fips, fips))

    def loadTracts(self):
        """Loads census tracts to PostGres
        Note: use cartographic boundary files, as "regular" tracts include coastline"""
        if 'tracts' in self.db.list_tables() and not self.forceUpdate:
            return

        tractFns = [ff for ff in os.listdir(paths['data']+'Tracts/') if ff.endswith('.shp') and ff.startswith('cb_2019')]
        self.db.execute('DROP TABLE IF EXISTS rawdata.tracts;')

        for ii, tractFn in enumerate(tractFns):
            # -a flag is omitted for first tract
            cmd = '''shp2pgsql -s 4326:4326 {} "{}Tracts/{}" rawdata.tracts | psql -d streetwidths -U adammb'''.format('-a'*(ii>0), paths['data'], tractFn)
            assert os.system(cmd) == 0
        self.db.execute('CREATE INDEX tracts_spat_idx ON rawdata.tracts USING gist (geom);')
        self.db.execute('CREATE INDEX tracts_spat_geogidx ON rawdata.tracts USING gist (geography(geom));') # to help distance matrix

        # drop water-only tracts
        self.db.execute('DELETE FROM rawdata.tracts WHERE aland=0;')

        # add census data
        df = self.getCensusDataFrame()
        df.drop(columns=['county','state','tract'], inplace=True)
        self.db.df2db(df, 'tracts_tmp', index=True)
        self.db.execute('''DROP TABLE IF EXISTS tracts2;
                           CREATE TABLE tracts2 AS
                           SELECT * from tracts LEFT JOIN tracts_tmp USING (geoid);''')
        self.db.execute('DROP TABLE tracts_tmp;')
        self.db.execute('DROP TABLE tracts;')
        self.db.execute('ALTER TABLE tracts2 RENAME TO tracts;')

        # create distance matrix for interpolation
        fipslist = ','.join(["'"+ff+"'" for ff in fipslookup.values()])
        cmd = """DROP TABLE IF EXISTS rawdata.tract_distances;
                CREATE TABLE rawdata.tract_distances AS
                SELECT t1.geoid geoid1, t2.geoid geoid2, ST_Distance(t1.geom::geography, t2.geom::geography) AS dist_m 
                FROM rawdata.tracts t1, rawdata.tracts t2 
                WHERE ST_DWithin(t1.geom::geography, t2.geom::geography, 5000) 
                      AND (t1.statefp || t1.countyfp IN ({fipslist})) AND (t2.statefp || t2.countyfp IN ({fipslist}));
                CREATE INDEX tract_distances_idx1 ON rawdata.tract_distances (geoid1);""".format(fipslist=fipslist)
        self.db.execute(cmd)

    def loadRoads(self):
        """Load OSM ways"""
        tables = self.db.list_tables()
        if all(['osm_'+fips in tables for fips in fipslookup.values()]) and not self.forceUpdate:
            return

        states = {'arizona':'04', 'california':'06','florida':'12', 'illinois':'17', 'massachusetts':'25', 'minnesota':'27',
                    'new-york':'36', 'tennessee': '47', 'texas':'48','washington':'53'}

        for state, fips in states.items():
            cmd = "java -Xmx5g -jar '{}osm2po-5/osm2po-core-5.2.43-signed.jar' tileSize=x cmd=c workDir='{}osm/' prefix='osm_{}' '{}Roads/{}-latest.osm.pbf'".format(paths['code'], paths['scratch'], fips, paths['data'], state)
            assert os.system(cmd)==0
            cmd = """psql -d streetwidths -h localhost -U adammb -q -f '{}osm/osm_{}_2po_4pgr.sql'""".format(paths['scratch'], fips)
            assert os.system(cmd)==0
            self.db.execute('CREATE INDEX osm_{}_2po_4pgr_spat_idx ON public.osm_{}_2po_4pgr USING gist (geom_way);'.format(fips, fips))

        shutil.rmtree(paths['scratch']+'osm')

        # Create county-specific roads file with projections
        for county, fips in fipslookup.items():
            srid = projections[county]
            self.db.execute('DROP TABLE IF EXISTS rawdata.osm_{};'.format(fips))
            self.db.execute('''CREATE TABLE rawdata.osm_{} AS
                               SELECT DISTINCT ON (w.id) w.* FROM public.osm_{}_2po_4pgr w, rawdata.tracts t
                               WHERE t.statefp = '{}' AND t.countyfp = '{}' 
                                     AND ST_Intersects(w.geom_way, t.geom);'''.format(fips, fips[:2], fips[:2], fips[2:]))
            self.db.execute("SELECT AddGeometryColumn('rawdata', 'osm_{}', 'geom', {},'LineString',2);".format(fips, srid))
            self.db.execute('UPDATE rawdata.osm_{} SET geom = ST_Transform(geom_way, {});'.format(fips, srid))
            self.db.execute('ALTER TABLE rawdata.osm_{} RENAME geom_way TO geom_wgs;'.format(fips))
            self.db.execute('CREATE INDEX osm_{}_source_idx ON rawdata.osm_{} (source);'.format(fips, fips))
            self.db.execute('CREATE INDEX osm_{}_target_idx ON rawdata.osm_{} (target);'.format(fips, fips))
            self.db.execute('CREATE INDEX osm_{}_spat_idx ON rawdata.osm_{} USING gist (geom);'.format(fips, fips))

        # drop original osm tables
        for state, fips in states.items():
            self.db.execute('DROP TABLE public.osm_{}_2po_4pgr;'.format(fips))

    def loadOther(self):
        """Loads urbanized area, raster of first built up year,
        and maybe other things that will get added later on"""
        if 'urbanized_area' in self.db.list_tables():
            return
        cmd = '''shp2pgsql -s 4326:4326 "{}Other/tl_2019_us_uac10.shp" rawdata.urbanized_area | psql -d streetwidths -U adammb'''.format(paths['data'])
        assert os.system(cmd) == 0
        self.db.execute('CREATE INDEX urbanized_area_spat_idx ON rawdata.urbanized_area USING gist (geom);')        

        cmd = '''DROP TABLE IF EXISTS main.tractareas;
                 CREATE TABLE main.tractareas AS
                    SELECT tracts.geoid, Substring(tracts.geoid, 1, 5) AS fips, Null::float AS urbanarea, 
                           SUM(ST_Area(ST_Intersection(tracts.geom, ST_Buffer(ua.geom,0))::geography)) / SUM(ST_Area(tracts.geom::geography)) AS urbanarea_frc
                    FROM rawdata.tracts, rawdata.urbanized_area ua
                    WHERE ST_Intersects(tracts.geom, ua.geom)
                    AND tracts.statefp || tracts.countyfp IN ({})
                    GROUP BY tracts.geoid;'''.format(','.join(["'"+vv+"'" for vv in fipslookup.values()]))
        self.db.execute(cmd)
        self.db.execute('''UPDATE main.tractareas ta 
                            SET urbanarea = t.aland * ta.urbanarea_frc / 1e6
                            FROM rawdata.tracts t WHERE t.geoid=ta.geoid;''')

        # load FBUY raster
        # reprojection process was trial and error
        # note: srid 97965 inserted manually using INSERT statement here: https://spatialreference.org/ref/sr-org/7965/
        cmd = 'raster2pgsql -d -C -M -I -t 100x100 {}FBUY/data/FBUY.tif rawdata.fbuytmp | psql -d streetwidths -U adammb'.format(paths['data'])
        assert os.system(cmd) == 0
        cmd = '''DROP TABLE IF EXISTS rawdata.fbuy;
                CREATE TABLE rawdata.fbuy AS
                SELECT rid, ST_SetSRID(rast, 97965) AS rast FROM rawdata.fbuytmp;
                DROP TABLE rawdata.fbuytmp;
                SELECT AddRasterConstraints('rawdata'::name, 'fbuy'::name, 'rast'::name);
                create index fbuy_sp_idx ON rawdata.fbuy USING gist((ST_ConvexHull(rast)));'''
        self.db.execute(cmd)

        # load land values from Davis 2019 or 2020
        df = self.getLandValueDf()
        self.db.df2db(df, 'landvalues', index=True)

    def getLandValueDf(self):
        """
        Create a dataframe with land values for each census tract, interpolating where necessary
        """
        outFn = paths['working'] + 'landvalues.pandas'
        if os.path.exists(outFn) and not self.forceUpdate:
            return pd.read_pickle(outFn)

        #landprice = pd.read_excel(paths['data']+'/Other/Davis 2019 Price of residential land.xlsx', sheet_name='Panel Census Tracts', skiprows=1)
        #landprice = landprice[landprice.Year==2018]
        inFn = paths['data']+'Other/Land-Prices_DLOS_2020_26Oct.xlsx'
        landprice = pd.read_excel(inFn, sheet_name='Cross-Section Census Tracts', skiprows=1)
        landprice['geoid'] = landprice['Census Tract'].astype(str).str.zfill(11)
        landprice['fips'] = landprice.geoid.str.slice(stop=5)
        landprice = landprice[landprice.fips.isin(fipslookup.values())]  # restrict to counties in our dataset

        # note: standardized is to 1/4 acre lots
        landprice.rename(columns={'Land Value\n(Per Acre, As-Is)':'landvalue_acre', "Land Value\n(1/4 Acre Lot, Standardized)":'landvalue_acre_std'}, inplace=True)
        landprice.landvalue_acre_std = landprice.landvalue_acre_std * 4  # convert to acre
        landprice.set_index('geoid', inplace=True)
        landprice = landprice[['landvalue_acre','landvalue_acre_std','fips']]

        # for interpolation: add lookup of nearby tracts (within 1km)
        tractLookUp = self.db.execfetchDf('SELECT * FROM rawdata.tract_distances;').set_index('geoid1')

        # restrict to tracts that actually have a landvalue
        tractLookUp = tractLookUp[tractLookUp.geoid2.isin(landprice.index.values)]

        # for interpolation: add lookup of nearby tracts (within 1km)
        interpGeoids = [gg for gg in tractLookUp.index.unique().values if gg not in landprice.index.values]
        interpPrices = {}
        for geoid in interpGeoids:
            near_geoids = tractLookUp.loc[geoid]
            if isinstance(near_geoids, pd.Series):  # with only 1 neighbor
                near_geoids = near_geoids.to_frame().T
            near_geoids = near_geoids.set_index('geoid2').join(landprice)
            assert len(near_geoids)>0
            # if any tracts are touching, just use those
            if any(near_geoids.dist_m==0):
                near_geoids = near_geoids[near_geoids.dist_m==0]  
            interpPrices[geoid] = near_geoids[['landvalue_acre','landvalue_acre_std']].mean()
        
        # merge back
        interpPrices = pd.DataFrame(interpPrices).T
        interpPrices.index.name = 'geoid'
        interpPrices['interpolated'] = True
        interpPrices['fips'] = interpPrices.index.str.slice(stop=5)
        landprice['interpolated'] = False
        landprice = pd.concat([landprice, interpPrices])

        # uprate prices to 2019 values based on the county-level change
        countyPrices = pd.read_excel(inFn, sheet_name='Panel Counties', skiprows=1)
        countyPrices['FIPS'] = countyPrices['FIPS'].astype(str).str.zfill(5)
        countyPrices.set_index('FIPS', inplace=True)
        countyPrices.rename(columns={'Land Value\n(Per Acre, As-Is)':'landvalue_acre', "Land Value\n(1/4 Acre Lot, Standardized)":'landvalue_acre_std'}, inplace=True)
        countyPrices = countyPrices[['Year','landvalue_acre','landvalue_acre_std']]
        countyPriceChange = countyPrices.loc[countyPrices.Year==2019] / countyPrices.loc[countyPrices.Year==2015]
        landprice = landprice.merge(countyPriceChange, left_on='fips',right_index=True, suffixes=('','_change'))
        for col in ['landvalue_acre','landvalue_acre_std']:
            landprice[col] = landprice[col]*landprice[col+'_change']
        
        landprice = landprice[['landvalue_acre','landvalue_acre_std','interpolated']]
        landprice.to_pickle(outFn)

        return landprice

    def getCensusDataFrame(self):
        """returns a df of the tables of interest for the counties of interest"""

        outFn = paths['working']+'census_compiled.pandas'
        if os.path.exists(outFn) and not self.forceUpdate:
            return pd.read_pickle(outFn) 

        statesToDo = np.unique([fp[:2] for fp in fipslookup.values()])
        state2letters = {'04':'az','06':'ca','12':'fl','17':'il','25':'ma','27':'mn','36':'ny','47':'tn','48':'tx','53':'wa',}
        path = paths['data']+'Census/'

        colsToUse = {'B01003_001':'totalpop', 'B11001_001':'totalhhs', 'B25001_001':'n_units',
                    'B25034_001':'n_units_check', 'B25034_002':'built_2014', 'B25034_003':'built_2010-13', 
                    'B25034_004':'built_2000-09', 'B25034_005':'built_1990-99', 'B25034_006':'built_1980-89',
                    'B25034_007':'built_1970-79', 'B25034_008':'built_1960-69', 'B25034_009':'built_1950-59',
                    'B25034_010':'built_1940-49', 'B25034_011':'built_pre-1940', 'B25035_001':'built_median'}
        geoColsToUse = ['LOGRECNO','COUNTY','STATE','TRACT','GEOID']
        geoFile, geoHeader = 'g20185XX.csv', '2018_SFGeoFileTemplate.xlsx'
        filesToDo = ['e20185XX0003000.txt','e20185XX0036000.txt','e20185XX0112000.txt', 'e20185XX0114000.txt']
        headers = ['seq3.xlsx','seq36.xlsx','seq112.xlsx','seq114.xlsx',]
        geoHeaders = pd.read_excel(path+'2018_5yr_Summary_FileTemplates/'+geoHeader,usecols=list(range(15))+[48,49])
        fileHeaders = [pd.read_excel(path+'2018_5yr_Summary_FileTemplates/'+hh) for hh in headers]

        statedfs = {}
        for statefips in statesToDo:
            state2letter = state2letters[statefips]
            counties = [cc[2:] for cc in fipslookup.values() if cc[:2]==statefips]
            geo = pd.read_csv(path+state2letter+'/'+geoFile.replace('XX',state2letter),header=None, usecols=list(range(15))+[48,49], encoding='latin-1')
            geo.columns=geoHeaders.columns

            # restrict to tracts in counties of interest
            geo = geo[geo.SUMLEVEL==140]
            geo = geo[geo.COUNTY.isin(counties)]
            geo = geo[geoColsToUse]

            # merge other files
            dfs = [pd.read_csv(path+state2letter+'/'+ff.replace('XX',state2letter), header=None) for ff in filesToDo]
            for df,header in zip(dfs, fileHeaders):
                df.columns = header.columns
                df.set_index('LOGRECNO', inplace=True)
            dfs = [df[[cc for cc in df.columns if cc in colsToUse]] for df in dfs]
            statedfs[statefips] = geo.set_index('LOGRECNO').join(dfs)   #[geoColsToUse+list(colsToUse.keys())]
            statedfs[statefips].rename(columns=colsToUse, inplace=True)

        # standardize geoid to match other datasets
        outDf = pd.concat(statedfs.values())
        outDf.columns = [cc.lower() for cc in outDf.columns]
        outDf['geoid'] = outDf.geoid.str.slice(start=7)
        outDf.set_index('geoid', inplace=True)
        outDf.county = outDf.county.astype(int).astype(str).str.zfill(3)
        outDf.state = outDf.state.astype(int).astype(str).str.zfill(2)
        outDf.tract = outDf.tract.astype(int).astype(str).str.zfill(6)

        outDf.loc[outDf.built_median.isin(['0','.']), 'built_median'] = np.nan
        outDf.built_median = outDf.built_median.astype(float)

        # check n_units
        assert all(outDf.n_units==outDf.n_units_check)
        outDf.drop(columns='n_units_check', inplace=True)

        # add land area and density
        landareas = pd.DataFrame(self.db.execfetch('SELECT geoid, aland/1e6 FROM rawdata.tracts;'), columns=['geoid','aland_sqkm']).set_index('geoid')
        outDf = outDf.join(landareas)
        outDf['units_sqkm'] = outDf.n_units.divide(outDf.aland_sqkm)
        outDf['pop_sqkm'] = outDf.totalpop.divide(outDf.aland_sqkm)

        outDf.to_pickle(outFn)
        return outDf

    def getHousingUnits(self):
        """Returns dataframe of urbanized housing units by county"""
        units = pd.read_csv(paths['data']+'Census/Housing_urbanrural_census2010.csv', skiprows=1, usecols=['id','Total','Total!!Urban!!Inside urbanized areas',])
        units.columns = ['id','total_units','urbanized_units']
        units['fips'] = units['id'].str.slice(start=9)
        units.set_index('fips', inplace=True)

        return units

    def loadAll(self):
        self.loadTracts()
        self.loadParcelsandRoads()
        self.loadRoads()
        self.loadOther()
        _ = self.getCensusDataFrame()
        _ = self.getLandValueDf()

class createStreetPolygons():
    """Creates polygons of streets based on the voids between parcels, using Voronoi polygons"""
    def __init__(self, forceUpdate=False):
        self.logger = logger('voronoi')
        self.db = tools.dbConnection(user='adammb', db='streetwidths', host='localhost', schema='main', requirePassword=False, verbose=True, logger=self.logger)
        self.db_rawdata = tools.dbConnection(user='adammb', db='streetwidths', host='localhost', schema='rawdata', requirePassword=False, verbose=True, logger=self.logger)
        self.forceUpdate = forceUpdate
        self.minEdgeLength = 15  # how long does an edge have to be to be included?

    def doCounty(self, fips):
        """Runs the analysis for one county"""

        if 'streetpoly_{}'.format(fips) in self.db.list_tables() and not self.forceUpdate:
            return

        print('Creating street polygons for {}'.format(fips))
        srid = projections[countylookup[fips]]

        # create temporary tracts in current SRID
        cmd = '''DROP TABLE IF EXISTS tmptracts_{fips};
                CREATE TABLE tmptracts_{fips} AS
                SELECT geoid, ST_Transform(geom, {srid}) AS geom
                FROM rawdata.tracts 
                WHERE countyfp='{countyfips}' AND statefp='{statefips}';
                CREATE INDEX tmptracts_{fips}_spidx ON main.tmptracts_{fips} USING GIST (geom);'''.format(fips=fips, srid=srid, statefips=fips[:2], countyfips=fips[2:])
        self.db.execute(cmd)

        # clip lines to parcel boundaries. 
        # This drops alleys, service roads, etc. that are mostly within the parcel polygon
        # buffer by -1 meter to account for inaccuracies (but also include parcels within 2.5 meters to avoid streets running along contigouous parcel boundaries)
        # some lines have discrete segments as a result, so we add the sid (segment id) unique field
        cmd = '''DROP TABLE IF EXISTS main.osm_clip_{fips};
            CREATE TABLE main.osm_clip_{fips} AS 
            SELECT id, osm_id, clazz, source, target, False as deadend, (ST_Dump(geom)).geom AS geom FROM
                (SELECT l.id, l.osm_id, l.clazz, l.source, l.target, 
                    CASE WHEN u.geom is Null OR ST_Area(u.geom)<2 THEN l.geom
                        ELSE ST_Difference(l.geom, ST_Buffer(u.geom,-1)) END AS geom
                FROM rawdata.osm_{fips} l 
                LEFT JOIN 
                    (SELECT l.id, ST_Union(p.geom) geom
                    FROM rawdata.osm_{fips} l, rawdata.parcels_{fips} p
                    WHERE ST_DWithin(l.geom, p.geom, 2.5)
                    GROUP BY l.id) u
                ON (u.id = l.id)) t1;
            ALTER TABLE main.osm_clip_{fips} ADD COLUMN sid SERIAL PRIMARY KEY;'''.format(fips=fips)
        self.db.execute(cmd)

        # identify deadends
        cmd = '''WITH nedges AS (
                SELECT node, COUNT(*) as nedges FROM
                    (SELECT source AS node FROM rawdata.osm_{fips}
                    UNION ALL 
                    SELECT target AS node FROM rawdata.osm_{fips}) t1
                GROUP BY node)
           UPDATE main.osm_clip_{fips} o
                SET deadend = (s.nedges = 1 OR t.nedges=1)
                FROM nedges s, nedges t
                WHERE o.source = s.node AND o.target = t.node;'''.format(fips=fips)
        self.db.execute(cmd)

        # convert lines into points in order to enable Voronoi polygon creation
        # ignore streets <= x meters
        cmd = '''DROP TABLE if exists main.streetpts_{fips};
            CREATE TABLE main.streetpts_{fips} AS
            WITH lines AS 
                   (SELECT sid, id, osm_id, clazz, deadend, geom FROM main.osm_clip_{fips} WHERE ST_Length(geom)>{elen}),
                ptfrcs AS
                   (SELECT sid, generate_series(5, ST_Length(geom)::int-5, 5)/ST_Length(geom) AS frc FROM lines)
            SELECT s.sid, s.id, s.osm_id, frc, s.clazz, s.deadend, 
                   Null::float as dist, Null::integer AS polyid, -- will be populated when calculating width
                   Null::smallint AS fbuy,
                   ST_LineInterpolatePoint(s.geom, frc) AS geom
            FROM  lines s, ptfrcs pf
            WHERE s.sid = pf.sid;'''.format(fips=fips, elen=self.minEdgeLength)
        self.db.execute(cmd)

        # add built-up year from raster
        cmd = '''UPDATE main.streetpts_{fips} t0
                SET fbuy = ST_Value(rast, ST_Transform(geom,97965)) 
                FROM rawdata.fbuy
                WHERE ST_Intersects(rast, ST_Transform(geom,97965));'''.format(fips=fips)
        self.db.execute(cmd)
        self.db.execute('UPDATE main.streetpts_{fips} SET fbuy = Null WHERE fbuy<2;'.format(fips=fips))

        # create voronoi polygons from points
        cmd = '''DROP TABLE IF EXISTS main.voronoi_{fips};
            CREATE TABLE main.voronoi_{fips} AS
            SELECT ST_MakeValid((st_dump(ST_VoronoiPolygons(geom))).geom) AS geom FROM 
                (SELECT ST_Collect(geom) as geom from main.streetpts_{fips}) p;
            ALTER TABLE main.voronoi_{fips} ADD COLUMN vid SERIAL PRIMARY KEY;
            CREATE INDEX voronoi_{fips}_spidx ON main.voronoi_{fips} USING GIST (geom);'''.format(fips=fips)
        self.db.execute(cmd)

        # dissolve voronoi polygons and add id of streets
        cmd = '''DROP TABLE IF EXISTS main.voronoi_diss_{fips};
            CREATE TABLE main.voronoi_diss_{fips} AS
            SELECT s.sid, s.id, s.osm_id, s.clazz, s.deadend, COUNT(*) as nvornoi, ST_Union(v.geom) AS geom FROM 
                main.voronoi_{fips} v,
                main.streetpts_{fips} s
            WHERE ST_Intersects(s.geom, v.geom)
            GROUP BY s.sid, s.id, s.osm_id, s.clazz, s.deadend;
            CREATE INDEX voronoi_diss_{fips}_spidx ON main.voronoi_diss_{fips} USING GIST (geom);'''.format(fips=fips)
        self.db.execute(cmd)

        # dissolve to OSM id (this avoids segmentation of polygons by lots of little side streets)
        # then dump in case polygons are not adjacent
        cmd = '''DROP TABLE IF EXISTS main.voronoi_diss2_{fips};
            CREATE TABLE main.voronoi_diss2_{fips} AS
            SELECT osm_id, clazz, bool_or(deadend) AS deadend, SUM(nvornoi) AS nvornoi, (ST_Dump(ST_MakeValid(ST_Union(geom)))).geom AS geom 
            FROM main.voronoi_diss_{fips} v
            GROUP BY osm_id, clazz;
            CREATE INDEX voronoi_diss2_{fips}_spidx ON main.voronoi_diss2_{fips} USING GIST (geom);'''.format(fips=fips)
        self.db.execute(cmd)  

        # clip to census tracts
        cmd = '''DROP TABLE IF EXISTS main.voronoi_clip_{fips};
            CREATE TABLE main.voronoi_clip_{fips} AS
            SELECT osm_id, clazz, deadend, nvornoi, ST_Intersection(v.geom, countygeom) AS geom 
            FROM 
                main.voronoi_diss2_{fips} v,
                (SELECT ST_Union(geom) AS countygeom FROM tmptracts_{fips}) t;
            CREATE INDEX voronoi_clip_{fips}_spidx ON main.voronoi_clip_{fips} USING GIST (geom);'''.format(fips=fips, statefips=fips[:2], countyfips=fips[2:])
        self.db.execute(cmd)    

        # add yearbuilt field if available
        if 0: # deprecated in favor if HISDAC
            cmd = '''ALTER TABLE main.voronoi_clip_{fips} ADD COLUMN yr_min float, ADD COLUMN yr_5pc float;'''.format(fips=fips)
            self.db.execute(cmd)        
            if 'yearbuilt' in self.db_rawdata.list_columns_in_table('parcels_{}'.format(fips)):
                cmd = '''UPDATE main.voronoi_clip_{fips} c1
                    SET yr_min = c2.yr_min, yr_5pc = c2.yr_5pc
                    FROM (SELECT v.osm_id, MIN(p.yearbuilt::integer) yr_min, percentile_disc(0.05) WITHIN GROUP (ORDER BY p.yearbuilt::integer) yr_5pc 
                            FROM rawdata.parcels_{fips} p, main.voronoi_clip_{fips} v
                            WHERE ST_Intersects(p.geom, v.geom) 
                            AND p.yearbuilt::integer IS NOT Null AND p.yearbuilt::integer>1800
                            GROUP BY v.osm_id) c2
                    WHERE c1.osm_id = c2.osm_id;'''.format(fips=fips)
                self.db.execute(cmd)        

        # remove part of polygons that overlaps with parcels to create ROW polygons
        # These will all be singlepart polygons, with a unique polyid field
        cmd = '''DROP TABLE IF EXISTS main.streetpoly_{fips};
            CREATE TABLE main.streetpoly_{fips} AS
            SELECT osm_id, clazz, deadend, nvornoi, (ST_Dump(geom)).geom AS geom 
            FROM (
                SELECT v.osm_id, v.clazz, v.deadend, v.nvornoi,  
                        CASE WHEN parcelgeom is Null THEN v.geom
                             WHEN {fips}=53033 AND v.osm_id=6522426 THEN ST_Difference(ST_buffer(v.geom,0.00001), parcelgeom) -- ugly kludge for one nasty parcel that throws an exception
                            ELSE ST_Difference(v.geom, parcelgeom) END AS geom
                FROM 
                    (SELECT osm_id, clazz, deadend, nvornoi, (ST_Dump(geom)).geom AS geom FROM main.voronoi_clip_{fips}) v
                    LEFT JOIN 
                    (SELECT v.osm_id, ST_Union(ST_MakeValid(p.geom)) AS parcelgeom 
                        FROM rawdata.parcels_{fips} p, main.voronoi_clip_{fips} v
                        WHERE ST_Intersects(p.geom, v.geom)
                        GROUP BY v.osm_id) u
                    USING (osm_id)
                WHERE ST_IsValid(v.geom)) t0;
            ALTER TABLE main.streetpoly_{fips} ADD COLUMN polyid SERIAL PRIMARY KEY;
            CREATE INDEX streetpoly_{fips}_spidx ON main.streetpoly_{fips} USING GIST (geom);'''.format(fips=fips)
        self.db.execute(cmd)

        # add estimate of average width
        # https://gis.stackexchange.com/questions/20279/calculating-average-width-of-polygon
        cmd = '''ALTER TABLE main.streetpoly_{fips} ADD COLUMN av_width float, ADD COLUMN av_width2 float, 
                       ADD COLUMN median_width float, ADD COLUMN length_m float, ADD COLUMN npts integer;
                UPDATE main.streetpoly_{fips}
                SET av_width = (ST_Perimeter(geom) - sqrt(POWER(ST_Perimeter(geom), 2) - 16 * ST_Area(geom)))/4
                WHERE POWER(ST_Perimeter(geom), 2) - 16 * ST_Area(geom)>0;'''.format(fips=fips)
        self.db.execute(cmd)

        # based on street centerline
        cmd = '''UPDATE main.streetpoly_{fips} t0
                SET av_width2 = CASE WHEN t1.length = 0 THEN Null 
                                ELSE t1.area / t1.length END 
                FROM
                    (SELECT p.polyid, ST_Area(p.geom) area,
                        SUM(ST_Length(ST_Intersection(l.geom, ST_Force2D(p.geom)))) length
                    FROM rawdata.osm_{fips} l, main.streetpoly_{fips} p
                    WHERE p.osm_id = l.osm_id
                    GROUP BY p.polyid) t1
                WHERE t0.polyid = t1.polyid;'''.format(fips=fips)
        self.db.execute(cmd)

        # median width - most robust
        # https://gis.stackexchange.com/questions/364173/calculate-median-width-of-road-polygons/364550#364550
        # assign this to point first, which allows for flexible aggregation to id or osm_id
        cmd = '''WITH frcs AS
                    (SELECT generate_series(0, 0.9, 0.2) AS frc1, generate_series(0.2, 1, 0.2) AS frc2),
            rings AS -- exterior ring of polygon, split into 5 substrings so that both sides of the street are captured separately
                (SELECT osm_id, polyid, 
                    CASE WHEN ST_NumInteriorRings(geom)=0 THEN ST_LineSubstring(ST_ExteriorRing(geom), frc1, frc2) -- second case is for circular roads
                    ELSE ST_LineSubstring(ST_Union(ST_ExteriorRing(geom),ST_InteriorRingN(geom, 1)), frc1, frc2) END AS geom 
                FROM frcs, main.streetpoly_{fips})
            UPDATE main.streetpts_{fips} t0
            SET dist = t3.avgdist, polyid = t3.polyid FROM
                (SELECT polyid, sid, frc, AVG(dist) avgdist -- average distance to the two closest segments of the polygon ring
                FROM
                    (SELECT ROW_NUMBER() OVER (PARTITION BY polyid, sid, frc ORDER BY dist) AS r, *
                    FROM 
                        (SELECT pt.osm_id, pt.polyid, pt.sid, pt.frc, ST_Distance(pt.geom, rings.geom) dist 
                            FROM rings,
                              (SELECT poly.polyid, pt0.osm_id, pt0.sid, pt0.frc, pt0.geom FROM main.streetpts_{fips} pt0,
                                 main.streetpoly_{fips} poly
                                 WHERE ST_Intersects(pt0.geom, poly.geom)) pt
                            WHERE rings.polyid = pt.polyid) t1) t2
                WHERE r<3 GROUP BY polyid, sid, frc) t3 -- limit to 2 closest segments
            WHERE t0.sid = t3.sid AND t0.frc = t3.frc;'''.format(fips=fips)
        self.db.execute(cmd)

        # set as zero when polygon intersects edge of county
        cmd = '''WITH b AS 
                    (SELECT ST_Transform(ST_ExteriorRing(ST_Union(geom)), {srid}) AS tractgeom
                     FROM rawdata.tracts WHERE statefp='{statefips}' AND countyfp='{countyfips}')
                UPDATE main.streetpts_{fips} pts
                SET dist = Null
                FROM  main.streetpoly_{fips} poly, b
                WHERE ST_Intersects(poly.geom, b.tractgeom)
                AND pts.polyid = poly.polyid;'''.format(fips=fips, srid=srid, statefips=fips[:2], countyfips=fips[2:])
        self.db.execute(cmd)

        # add median to polygons
        cmd = '''UPDATE main.streetpoly_{fips} t0           
                 SET median_width = t1.width, npts=t1.npts  FROM
                    (SELECT polyid, COUNT(*) as npts, percentile_cont(0.5) within group (order by dist)*2 AS width -- median distance * 2
                        FROM main.streetpts_{fips} 
                        GROUP BY polyid) t1
                 WHERE t0.polyid = t1.polyid;'''.format(fips=fips)
        self.db.execute(cmd)

        # create osm_id level table with median width and length
        cmd = '''DROP TABLE IF EXISTS main.osm_metrics_{fips};
                 CREATE TABLE main.osm_metrics_{fips} AS
                 SELECT osm_id, bool_or(deadend) AS deadend, max(clazz) as clazz, --clazz are all the same
                        COUNT(*) as npts, percentile_cont(0.5) within group (order by dist)*2 AS median_width, -- median distance * 2
                        Null::float as length, Null::float as netlength, Null::float AS area, Null::float AS perimeter,
                        Null::float as intersection_width, Null::bool as urbanized, 
                        percentile_disc(0.5) WITHIN GROUP (ORDER BY fbuy) AS fbuy,
                        FROM main.streetpts_{fips} 
                        GROUP BY osm_id;'''.format(fips=fips)
        self.db.execute(cmd)

        # add intersection widths 
        cmd = '''WITH w AS
                    (SELECT o.osm_id, o.source, o.target, m.median_width
                        FROM rawdata.osm_{fips} o
                        LEFT JOIN main.osm_metrics_{fips} m USING (osm_id))
            UPDATE main.osm_metrics_{fips} metrics
            SET intersection_width = t2.intersection_width
            FROM (
                SELECT osm_id, SUM(width) AS intersection_width -- sum over each end of each edge
                FROM (
                    SELECT id, osm_id, sourcenode, AVG(width) AS width 	-- width of half of each intersection
                        FROM (
                        SELECT o.id AS id, o.osm_id AS osm_id, median_width/2. AS width, TRUE as sourcenode
                            FROM rawdata.osm_{fips} o, w 
                            WHERE (o.source=w.source OR o.source=w.target) AND o.osm_id!=w.osm_id
                        UNION ALL
                        SELECT o.id AS id, o.osm_id AS osm_id, median_width/2. AS width, FALSE as sourcenode
                            FROM rawdata.osm_{fips} o, w 
                            WHERE (o.target=w.source OR o.target=w.target) AND o.osm_id!=w.osm_id) t0
                    GROUP BY id, osm_id, sourcenode) t1
                GROUP BY osm_id) t2
            WHERE metrics.osm_id = t2.osm_id;'''.format(fips=fips)
        self.db.execute(cmd)

        # add street length to polygons
        cmd = '''UPDATE main.osm_metrics_{fips} metrics           
                 SET length = t1.length FROM
                    (SELECT osm_id, SUM(ST_Length(geom)) AS length
                      FROM main.osm_clip_{fips}
                      GROUP BY osm_id) t1
                WHERE metrics.osm_id = t1.osm_id;'''.format(fips=fips)
        self.db.execute(cmd) 

        cmd = '''UPDATE main.osm_metrics_{fips}           
                 SET netlength = GREATEST(0, length - COALESCE(intersection_width,0));'''.format(fips=fips)
        self.db.execute(cmd)
    
        cmd = '''UPDATE main.osm_metrics_{fips} metrics
               SET area = poly.area, perimeter = poly.perimeter FROM
                  (SELECT osm_id, SUM(ST_Area(geom)) AS area, SUM(ST_Perimeter(geom)) as perimeter
                       FROM main.streetpoly_{fips}
                       GROUP BY osm_id) poly
               WHERE metrics.osm_id = poly.osm_id'''.format(fips=fips)
        self.db.execute(cmd)

        # add urbanized flag. buffer is to avoid topological issues
        cmd = '''ALTER TABLE main.streetpoly_{fips} ADD COLUMN urbanized bool;
                 WITH ua AS (
                     SELECT ST_Union(st_buffer(ST_Transform(u.geom, {srid}),0.1)) AS geom 
                        FROM rawdata.urbanized_area u, 
                             (SELECT ST_Collect(geom) AS tractgeom
                                 FROM rawdata.tracts WHERE statefp='{statefips}' AND countyfp='{countyfips}') t
                        WHERE ST_Intersects(u.geom, t.tractgeom))
                 UPDATE main.streetpoly_{fips} s
                 SET urbanized = ST_Intersects(s.geom, ua.geom)
                 FROM ua;'''.format(fips=fips, srid=srid, statefips=fips[:2], countyfips=fips[2:])
        self.db.execute(cmd)

        cmd = '''UPDATE main.osm_metrics_{fips} s
                 SET urbanized = t.urbanized -- , yr_min=t.yr_min, yr_5pc=t.yr_5pc
                 FROM
                    (SELECT osm_id, bool_or(urbanized) AS urbanized 
                        -- , AVG(yr_min) AS yr_min, AVG(yr_5pc) AS yr_5pc -- (deprecated) both will be the same within osm_id
                    FROM main.streetpoly_{fips} 
                    GROUP BY osm_id) t
                WHERE s.osm_id = t.osm_id;'''.format(fips=fips)
        self.db.execute(cmd)

        # lookup to census tract
        cmd = '''DROP TABLE IF EXISTS main.streettotract_{fips};
            CREATE TABLE main.streettotract_{fips} AS
            SELECT s.osm_id, t.geoid FROM 
                main.osm_clip_{fips} s,
                main.tmptracts_{fips} t
            WHERE ST_Intersects(s.geom, t.geom);'''.format(fips=fips, statefips=fips[:2], countyfips=fips[2:])
        self.db.execute(cmd)

        # add value measures
        cmd = '''DROP TABLE IF EXISTS tmp_osm_metrics_{fips};
            CREATE TABLE tmp_osm_metrics_{fips} AS 
            SELECT *, area / {sqm_per_acre} * landvalue_acre / 1e6 AS streetvalue_m,
                    area / {sqm_per_acre} * landvalue_acre_std / 1e6 AS streetvalue_m_std,
                    clazz=41 AS is_residential,
                    median_width * {m2ft} AS width_ft,
                        -- drop "irregular" streets where "regularized" area is less than 1/4 of actual area, or perimeter is long. Mostly, these have orphaned bits of land that are attached
                    (area/(length*median_width)>5) OR (perimeter/(length+median_width)>8) AS ignore
            FROM main.osm_metrics_{fips}
            LEFT JOIN 
            (
            SELECT osm_id, AVG(landvalue_acre) AS landvalue_acre, AVG(landvalue_acre_std) AS landvalue_acre_std,
                AVG(interpolated::int)::real AS interpolated, AVG(built_median) AS built_median,
                AVG(units_sqkm) AS units_sqkm, AVG(pop_sqkm) AS pop_sqkm 
            FROM main.streettotract_{fips} s2t
            LEFT JOIN rawdata.landvalues lv USING(geoid)
            LEFT JOIN rawdata.tracts t USING (geoid)
            GROUP BY osm_id) t0
            USING (osm_id);'''.format(fips=fips, sqm_per_acre=sqm_per_acre, m2ft=m2ft)
        self.db.execute(cmd)
        self.db.execute('ALTER TABLE tmp_osm_metrics_{fips} ADD COLUMN value_std_per_25ft float;'.format(fips=fips))
        self.db.execute('UPDATE tmp_osm_metrics_{fips} SET value_std_per_25ft = streetvalue_m_std*1e6/length/{m2ft}*25;'.format(fips=fips, m2ft=m2ft))
        self.db.execute('DROP TABLE osm_metrics_{fips};'.format(fips=fips))
        self.db.execute('ALTER TABLE tmp_osm_metrics_{fips} RENAME TO osm_metrics_{fips};'.format(fips=fips))

        # clean up
        for table in ['osm_clip','voronoi','voronoi_diss','voronoi_diss2','voronoi_clip','tmptracts']:  # streetpts
            self.db.execute('DROP TABLE main.{table}_{fips};'''.format(table=table, fips=fips))

    def doAllCounties(self):
        """Wrapper for doCounty"""
        for fips in fipslookup.values():
            self.doCounty(fips)

class censusGetter():
    """NOT USED"""
    def __init__(self):
        self.apiKey = 'fillinhere'

    def getACSTable(self, table, statefips, countyfips):
        """ DO NOT USE - not as much detail available as with web interface or FTP"""
        import requests
        base = 'https://api.census.gov/data/'
        year = '2018'
        apiDict = {'year':year, 'tablename':table, 'statefips':statefips, 'countyfips':countyfips, 'apiKey':self.apiKey}
        apiCall = base+'{year}/acs/acs5/subject?get={tablename}&for=tract:*&in=state:{statefips}%20county:{countyfips}&key={apiKey}'.format(**apiDict)
        response = requests.get(apiCall).json()
        df = pd.DataFrame(response[1:], columns=response[0])
        return df

    def downloadAllTables(self):
        tablesToGet = {'S0101_C01_001E':'totalpop', 'S1101_C01_001E':'totalhhs', 'S2501_C01_001E':'units',
            'S2504_C01_009E':'built2014-', 'S2504_C01_010E':'built2010-13', 'S2504_C01_011E':'built2000-09',
            'S2504_C01_012E':'built1980-99', 'S2504_C01_013E':'built1960-79', 
            'S2504_C01_014E':'built1940-59', 'S2504_C01_015E':'builtpre1940',}
        allTables = ','.join(tablesToGet.keys())
        df = pd.concat([self.getACSTable(allTables, fips[:2], fips[2:]) for fips in fipslookup.values()])
        df.rename(columns=tablesToGet, inplace=True)

class analysis():
    """Produces charts, etc. for paper"""

    def __init__(self, cutoff=25, bw=False):
        """
        Cutoff is minimum length of residential streets (25m)
        """
        self.db = tools.dbConnection(user='adammb', db='streetwidths', host='localhost', schema='main', requirePassword=False)
        self.bw = bw
        self.streetDf = None # populated in self.getDf
        self.tractareas = self.db.execfetchDf('SELECT * from main.tractareas;')
        self.countyareas = self.tractareas.groupby('fips').urbanarea.sum() 
        self.countyareas['all'] = self.countyareas.sum()
        self.cutoff=cutoff
        # from osm2po config file
        self.clazzDict = {11: 'motorway', 12: 'motorway_link', 13:'trunk', 14:'trunk_link', 
                        15:'primary', 16:'primary_link', 21:'secondary', 22:'secondary_link', 
                        31:'tertiary', 32:'tertiary_link', 41:'residential', 42:'road', 
                        43:'unclassified', 51:'service', 63:'living_street'}

    def getDf(self, urbanized=True, saveDf=False):
        """Gets dataframe of streets with widths and tract id from postgres\
            urbanized=True restrict to streets that intersect census-defined urbanized area
            
        Note: returned df is saved (for convenience) in Working/streets2.pandas"""

        if self.streetDf is not None:
            return self.streetDf

        print ('Compiling street dataframe')
        dfs, streets, streettotracts = {}, {}, {}
        for fips in fipslookup.values():
            streets[fips] = self.db.execfetchDf('SELECT {fips} as fips, osm_id, clazz, deadend, npts, median_width, length, netlength, area, perimeter, intersection_width, urbanized, fbuy, streetvalue_m, streetvalue_m_std, value_std_per_25ft, landvalue_acre, landvalue_acre_std, pop_sqkm, is_residential, width_ft, ignore FROM main.osm_metrics_{fips};'''.format(fips=fips))
            assert streets[fips]['osm_id'].is_unique
        
        streetDf = pd.concat(streets.values())
        del streettotracts, streets, dfs

        if 0: # this is now done in createStreetPolygons() 
            streettotract = pd.concat(streettotracts.values())
            streettotract['fips'] = streettotract.fips.astype(str).str.zfill(5)

            # now we can calculate value of streets! also add in covariates from tracts
            dl = dataLoader()
            landvalue = dl.getLandValueDf()
            tractdata = dl.getCensusDataFrame()

            streettotract = streettotract.set_index('geoid').join(landvalue).join(tractdata).reset_index()
            streettotract['interpolated'] = streettotract.interpolated.astype(float)  # 

            # issue is that some streets abut multiple tracts, so do a groupby and then join
            joinDf = streettotract.groupby('osm_id')[['landvalue_acre','landvalue_acre_std','interpolated','built_median','units_sqkm','pop_sqkm']].mean()

            streetDf = streetDf.set_index('osm_id').join(joinDf)

            streetDf['streetvalue_m'] = streetDf.area / sqm_per_acre * streetDf.landvalue_acre / 1e6
            streetDf['streetvalue_m_std'] = streetDf.area / sqm_per_acre * streetDf.landvalue_acre_std / 1e6
            streetDf['is_residential'] = streetDf.clazz==41
            streetDf['width_ft'] = streetDf.median_width * m2ft

        streetDf['fips'] = streetDf.fips.astype(str).str.zfill(5)
        streetDf.set_index('osm_id', inplace=True)

        if urbanized:
            streetDf = streetDf[streetDf.urbanized==True]

        # drop "irregular" streets where "regularized" area is less than 1/4 of actual area. Mostly, these have orphaned bits of land that are attached
        if 0: # now done in postgres
            len1 = len(streetDf)
            streetDf = streetDf[(streetDf.area/(streetDf.length*streetDf.median_width))<=5]
            len2 = len(streetDf)
            print('Dropped {} ({:.1f}%) irregular streets where area is too large'.format((len1-len2),(len1-len2)/len1*100 ))

            streetDf = streetDf[(streetDf.perimeter/(streetDf.length+streetDf.median_width))<=8]
            len3 = len(streetDf)
            print('Dropped {} ({:.1f}% of original sample) irregular streets where perimeter is too long'.format((len2-len3),(len2-len3)/len1*100 ))
        else:
            print('Dropping {} ({:.1f}%) irregular streets where area is too large'.format(streetDf.ignore.sum(),streetDf.ignore.sum()/len(streetDf)*100 ))
            print('({:.1f}%) of residential streets'.format(streetDf.loc[streetDf.is_residential,'ignore'].sum()/streetDf.is_residential.sum()*100 ))
            streetDf = streetDf[streetDf.ignore==False]

        self.streetDf = streetDf
        self.streetDf.urbanized = self.streetDf.urbanized.astype(bool)  # allows saving in Stata
        # also populate dataframe of urban residential streets  
        self.dfr = streetDf[(streetDf.urbanized==True) & (streetDf.clazz==41)]  # residential

        if saveDf:
            assert self.streetDf.ignore.sum()==0 and self.dfr.ignore.sum()==0 # we dropped all of them
            self.streetDf.drop(columns='ignore', inplace=True) # to enable Stata export
            self.dfr.drop(columns='ignore', inplace=True) # to enable Stata export
            self.streetDf.to_pickle(paths['working']+'streets_all.pandas')
            self.streetDf.to_stata(paths['working']+'streets_all.dta')
            self.dfr.to_pickle(paths['working']+'streets_residential.pandas')
            self.dfr.to_stata(paths['working']+'streets_residential.dta')

        return self.streetDf

    def widthHistograms(self):
        """
        Table 2 (street widths and area) and Fig 1 (histograms)
        """

        df = self.getDf().copy()
        assert df.urbanized.all()
        df['widthxlength'] = df.width_ft * df.netlength # for weighted sume
        aggDfs = {}
        for ii, groupCols in enumerate([['is_residential','fips'],['is_residential']]):
            aggDfs[ii] = df.groupby(groupCols).width_ft.agg(['median','mean'])
            aggDfs[ii].columns=['Median width (ft)','Mean width (ft)']
            aggDfs[ii]['Mean width (ft) (weighted)'] = df.groupby(groupCols).widthxlength.sum() / df.groupby(groupCols).netlength.sum()
            aggDfs[ii]['Area (km2)'] = df.groupby(groupCols).area.sum() / 1e6
            aggDfs[ii]['Length (km)'] = df.groupby(groupCols).netlength.sum() / 1e3
        aggDfs[1]['fips'] = 'all'
        aggDfs[1].set_index('fips', append=True, inplace=True)
        aggDf = pd.concat([aggDfs[0], aggDfs[1]])
        aggDf['Area (percent)'] = aggDf.apply(lambda x: x['Area (km2)'] / self.countyareas[x.name[1]], axis=1) *100

        tableDf = aggDf.loc[True][['Median width (ft)','Mean width (ft)','Mean width (ft) (weighted)','Length (km)','Area (km2)','Area (percent)']].round(1).astype(str)
        tableDf['All streets area (percent)'] = aggDf.groupby(level=1)['Area (percent)'].sum().round(1).astype(str)
        tableDf.reset_index(inplace=True)
        tableDf['County'] = tableDf.fips.map(countylookup_formatted)
        print(tableDf)
        tableDf.to_csv(paths['graphics']+'Table1_dimensions.csv', index=False)

        #### Fig 3  Distributions of width
        dfr = self.dfr.copy()
        dfr = dfr[(dfr.length>=self.cutoff) & (dfr.width_ft<=75)] # exclude from histogram streets that are short or have extreme widths
        assert dfr.urbanized.all()
        dfr.width_ft = dfr.width_ft.round()

        def plotHist(fig, axes, countiesToPlot, fname):
            # in function because we do a separate version for the appendix
            barcolor = '0.5'
            labelpos = (0.05, 0.85) if len(countiesToPlot)==4 else (0.05, 0.65)
            for ii, county in enumerate(countiesToPlot):
                if county=='all':
                    dfr.width_ft.plot(kind='hist',bins=range(15,76), ax=axes[0,0], align='left', color=barcolor)
                    axes[0,0].set_title('All counties', position=labelpos, ha='left')
                else:
                    row, col = math.floor(ii/2), ii%2
                    dfr[dfr.fips==county].width_ft.plot(kind='hist',bins=range(15,76), ax=axes[row,col], align='left', color=barcolor)
                    axes[row, col].set_title(countylookup_formatted[county], position=labelpos, ha='left', fontsize=10)

            for row in range(axes.shape[0]):
                axes[row,1].set_ylabel('')
                # for appendix plot, only have values on bottom row
                if len(countiesToPlot)>4 and row!=(axes.shape[0]-1):
                    axes[row,0].set_xticklabels([])
                    axes[row,1].set_xticklabels([])

            for col in [0,1]:
                axes[-1, col].set_xlabel('Width (feet)')

            for ax in axes.flatten():
                ax.set_xlim(15,76)
                ax.set_xticks(range(20,76,10))
                ax.set_xticks(range(15,76,1), minor=True)
                ax.set_yticklabels([])
            fig.tight_layout()
            fig.savefig(paths['graphics']+fname)

        fig, axes = plt.subplots(2,2, figsize = figsizes['full'])
        main_counties = ['all','04013','12086','06071',] # ['all','12086','06071','48439']  # counties to plot in main figure
        plotHist(fig, axes, main_counties, 'Fig3_width_distributions.pdf')
        fig, axes = plt.subplots(10,2, figsize = (7,10))
        plotHist(fig, axes, fipslookup.values(), 'FigA4_width_distributions.pdf')

        # get back the original with the areas of short and overly-wide streets
        dfr = self.dfr.set_index('fips').copy()   
        assert dfr.urbanized.all()
        dfr.width_ft = dfr.width_ft.round()

        print('Note: all figures are weighted by length')
        print('{:.2f} percent of streets are 50 or 60 ft wide'.format(dfr[dfr.width_ft.isin([50,60])].length.sum() / dfr.length.sum()*100))
        tmpDf = dfr.loc['48439']
        print('{:.2f} percent of Tarrant streets are 50 ft wide'.format((tmpDf[tmpDf.width_ft==50]).length.sum() / tmpDf.length.sum()*100))
        tmpDf = dfr.loc['04013']
        print('{:.2f} percent of Maricopa streets are 50 ft wide'.format((tmpDf[tmpDf.width_ft==50]).length.sum() / tmpDf.length.sum()*100))
        tmpDf = dfr.loc['17031']
        print('{:.2f} percent of Cook streets are 66 ft wide'.format((tmpDf[tmpDf.width_ft==66]).length.sum() / tmpDf.length.sum()*100))
        for yr in [1960,1970,1980,1990]:
            tmpDf = dfr.loc['12086'][dfr.loc['12086','fbuy']>=yr]
            print('{:.2f} percent of post-{} streets in Miami Dade are 50 ft wide'.format((tmpDf[tmpDf.width_ft==50]).length.sum() / tmpDf.length.sum()*100, yr))

    def timeTrends_old(self):
        """Dating measures using parcel and census year built"""
        _ = self.getDf()
        dfr = self.dfr.copy()

        def rollmean(adf, minyr, meanCol='width_ft', wtCol='netlength'):
            # 5-year weighted rolling mean. Hard to do this with pd.rolling
            adf = adf[adf.index.notnull()].copy()
            adf.sort_index(inplace=True)
            result = pd.Series({yr: (adf.loc[yr-5:yr+5, meanCol]*adf.loc[yr-5:yr+5, wtCol]).sum()/(adf.loc[yr-5:yr+5, wtCol]).sum() for yr in range(minyr,2020) if len(adf.loc[yr-5:yr+5])>=100})
            if minyr==1900: result[1900] =  (adf.loc[:1900, meanCol]*adf.loc[:1900, wtCol]).sum()/(adf.loc[:1900, wtCol]).sum()
            return result.sort_index()

        yrcounts = dfr.groupby('fips').yr_5pc.count()
        with_yrbuilt = yrcounts[yrcounts>1].index.values
        fig, ax = plt.subplots(figsize=figsizes['full'])
        rollmean(dfr.set_index('yr_5pc'), 1900).plot(ax=ax, label='All counties (parcel-based date)')
        rollmean(dfr.set_index('built_median'), 1940).plot(ax=ax, label='All counties (census-based date)')
        for fips in fipslookup.values():
            if fips in with_yrbuilt:
                rollmean(dfr[dfr.fips==fips].set_index('yr_5pc'), 1900).plot(ax=ax, color='0.6', lw=0.5, label='_nolegend_')
            elif fips in ['06085','53033']: # include in legend
                rollmean(dfr[dfr.fips==fips].set_index('built_median'), 1940).plot(ax=ax, lw=1,label=countylookup_formatted[fips]) 
            else:
                rollmean(dfr[dfr.fips==fips].set_index('built_median'), 1940).plot(ax=ax, color='0.6', lw=0.5,label='_nolegend_' ) 
        ax.set_ylim(42,64.5)
        ax.set_xlim(1900,2020)
        ax.set_xticklabels(['<1900']+[str(yr) for yr in range(1920,2021,20)], fontsize=12)
        ax.legend()
        ax.set_ylabel('Width (weighted mean, ft)', fontsize=12)
        fig.tight_layout()
        fig.savefig(paths['graphics']+'Fig2_trends.pdf')

        # what about standard deviation over time?
        def rollstd(adf, minyr, col='width_ft'):
            # 5-year weighted rolling mean. Hard to do this with pd.rolling
            adf = adf[adf.index.notnull()].copy()
            adf.sort_index(inplace=True)
            result = pd.Series({yr: adf.loc[yr-5:yr+5, col].std() for yr in range(minyr,2020) if len(adf.loc[yr-5:yr+5])>=150})
            if minyr==1900: result[1900] =  adf.loc[:1900, col].std()
            return result.sort_index()

        fig, ax = plt.subplots(figsize=figsizes['full'])
        rollstd(dfr.set_index('yr_5pc'), 1900).plot(ax=ax, label='All counties (parcel-based date)')
        rollstd(dfr.set_index('built_median'), 1940).plot(ax=ax, label='All counties (census-based date)')
        for fips in fipslookup.values():
            if fips in with_yrbuilt:
                rollstd(dfr[dfr.fips==fips].set_index('yr_5pc'), 1900).plot(ax=ax, color='0.6', lw=0.5, label='_nolegend_')
            elif fips in ['06085','53033']: # include in legend
                rollstd(dfr[dfr.fips==fips].set_index('built_median'), 1940).plot(ax=ax, lw=1,label=countylookup_formatted[fips]) 
            else:
                rollstd(dfr[dfr.fips==fips].set_index('built_median'), 1940).plot(ax=ax, color='0.6', lw=0.5,label='_nolegend_' ) 
        ax.set_xlim(1900,2020)
        ax.set_xticklabels(['<1900']+[str(yr) for yr in range(1920,2021,20)], fontsize=12)
        ax.legend()
        ax.set_ylabel('Width (weighted mean, ft)', fontsize=12)
        fig.tight_layout()
        fig.savefig(paths['graphics']+'Fig_trends_stddev.pdf')

    def timeTrends(self):
        """Dating measures using FBUY field from Leyk et al 2020"""
        _ = self.getDf()
        dfr = self.dfr.copy()
        if self.bw:
            colors = ['k']*3+['0.5']*3
            lstyles = ['dotted','dashed','dashdot']*2
        else:
            colors = c5s
            lstyles = ['solid']*6
        specialfips = {'06001':0,'06085':1,'12086':2,'48439':3,'53033':4,'47157':5}

        # consolidate pre-1850 development
        dfr = dfr[pd.notnull(dfr.fbuy)]
        dfr.loc[dfr.fbuy<1850, 'fbuy'] = 1850 
        dfr.set_index('fbuy', inplace=True)
        dfr.index.name = '' # avoids plotting xlabel

        fig, ax = plt.subplots(figsize=figsizes['full'])
        plotDf = (dfr.width_ft*dfr.netlength).groupby(level=0).sum()/dfr.groupby(level=0).netlength.sum()
        plotDf.plot(ax=ax, color='k', lw=3, label='All counties')

        for fips in fipslookup.values():
            mask = dfr.fips==fips
            plotDf = (dfr[mask].width_ft*dfr[mask].netlength).groupby(level=0).sum()/dfr[mask].groupby(level=0).netlength.sum()
            # limit to years with > n streets
            yrcounts = dfr[mask].groupby(level=0).size()
            plotDf = plotDf[yrcounts>=20]
            if 2015 in plotDf: print(fips, plotDf.loc[2015])
            if fips in specialfips: # include in legend
                plotDf.plot(ax=ax, lw=1.2,label=countylookup_formatted[fips], color=colors[specialfips[fips]], linestyle=lstyles[specialfips[fips]],) 
            else:
                plotDf.plot(ax=ax, color='0.6', lw=0.5, label='_nolegend_' ) 
        ax.set_ylim(33,64.5)
        ax.set_xticks(range(1850,2015,25))
        ax.set_xticklabels(['<1850']+[str(yy) for yy in range(1875,2015,25)], fontsize=12)
        ax.set_xlim(1850,2015)
        #ax.spines['right'].set_visible(False)
        #ax.spines['top'].set_visible(False) 
        ax.legend(ncol=2, loc=8)
        ax.set_ylabel('Width (weighted mean, feet)', fontsize=12)
        fig.tight_layout()
        fig.savefig(paths['graphics']+'Fig4_trends'+'_bw'*self.bw+'.pdf')

        # STANDARD DEVIATION OF WIDTH
        fig, ax = plt.subplots(figsize=figsizes['full'])
        dfr.groupby(level=0).median_width.std().plot(ax=ax, label='All counties')
        for fips in fipslookup.values():
            plotDf = dfr[dfr.fips==fips].median_width
            yrcounts = dfr[dfr.fips==fips].groupby(level=0).size()
            plotDf = plotDf[yrcounts>=20]
            if fips in ['06085','53033']: # include in legend
                plotDf.groupby(level=0).std().plot(ax=ax, lw=1,label=countylookup_formatted[fips]) 
            else:
                plotDf.groupby(level=0).std().plot(ax=ax, color='0.6', lw=0.5,label='_nolegend_' ) 
        ax.set_xticks(range(1850,2015,25))
        ax.set_xticklabels(['<1850']+[str(yy) for yy in range(1875,2015,25)], fontsize=12)
        ax.legend()
        ax.set_ylabel('Standard deviation of width', fontsize=12)
        fig.tight_layout()
        fig.savefig(paths['graphics']+'Fig_trends_stddev.pdf')

    def streetValue(self):
        """calculate value of rights of way"""
        df = self.getDf().copy()
        units = dataLoader().getHousingUnits()
        colors = ['0.1','0.7'] if self.bw else c2s
        
        aggDf = (df.groupby(['fips','is_residential']).streetvalue_m_std.sum()).to_frame()
        aggDf = aggDf.join(units['urbanized_units'])

        # add totals and county name
        aggDf2 = aggDf.groupby(level=1).sum().reset_index()
        aggDf2['county'] = 'All counties'
        aggDf.reset_index(inplace=True)
        aggDf['county'] = aggDf.fips.map(countylookup_formatted)
        aggDf = pd.concat([aggDf, aggDf2])

        # per unit values (thousand $)
        aggDf['streetvalue_1k_perunit'] = aggDf.streetvalue_m_std*1e3/aggDf.urbanized_units
        aggDf['is_residential'] = aggDf.is_residential.apply(lambda x: 'Residential streets' if x else 'Other streets')
        aggDf['streetvalue_bn_std'] = aggDf.streetvalue_m_std/1e3

        # reshape for plotting
        plotDf1 = aggDf[pd.notnull(aggDf.fips)].set_index(['county','is_residential']).streetvalue_bn_std.unstack()[['Residential streets','Other streets']]
        plotDf2 = aggDf.set_index(['county','is_residential']).streetvalue_1k_perunit.unstack()[['Residential streets','Other streets']]
        plotDf2=pd.concat([plotDf2.iloc[1:2], plotDf2.iloc[0:1], plotDf2.iloc[2:]])   # kludge for sort order
        fig, axes = plt.subplots(1,2, figsize=(7,5))  # 0.5 inch wider than "full"
        plotDf1.plot(kind='bar', stacked=True, ax=axes[0], legend=None, color=colors)
        axes[0].set_ylabel('Land value (billion $)', fontsize=11)
        axes[0].set_xlabel('') 
        axes[0].xaxis.set_tick_params(labelsize=9)

        plotDf2.plot(kind='bar', stacked=True, ax=axes[1], color=colors)
        axes[1].set_ylabel("Land value per housing unit (thousand $)", fontsize=11, labelpad=15, rotation=270)
        axes[1].set_xlabel('')
        axes[1].xaxis.set_tick_params(labelsize=9)
        axes[1].yaxis.tick_right()
        axes[1].yaxis.set_label_position("right")
        # flip legend
        handles, labels = axes[1].get_legend_handles_labels()
        axes[1].legend(reversed(handles), reversed(labels), title=None, loc='upper left', frameon=False)

        fig.tight_layout()

        fig.savefig(paths['graphics']+'Fig5_value_histogram'+'_bw'*self.bw+'.pdf')

        # save for Tableau chart
        aggDf.to_csv(paths['website'] + 'Data/street_values.tsv', sep='\t', index=False)

        # statistics for text
        print('Total value of streets: ${:.1f} billion'.format(plotDf1.sum().sum()))
        print('Total value of residential streets: ${:.1f} billion'.format(plotDf1['Residential streets'].sum()))
        for fips in ['48029','48439','47157','36047','36081','06037','06001','06085']:
            cname = countylookup_formatted[fips]
            print('Value of residential streets per unit in {}: ${:.0f},000'.format(cname, plotDf2.loc[cname, 'Residential streets']))

        # counterfactuals with reducing width
        # complication here is in allocating intersections, irregularities, etc. 
        # We do this by calculating % reduction in the "regular" area, and applying that proportionately

        dfr = self.dfr.set_index('fips').copy()   # so that we include the areas of short and long streets
        assert dfr.urbanized.all()

        valueRatio = 1 - dfr[pd.isnull(dfr.median_width)].streetvalue_m_std.sum() / dfr.streetvalue_m_std.sum()   # adjustment factor to take into account streets with no width
        maxW = 70
        widths = list(np.arange(0,maxW+1, 0.2))
        stvalues, stvalues_ph = {}, {}
        for fips in countylookup:
            # calculate value based on (constrained width / width) * value
            stvalues[fips] = [(np.minimum(dfr.loc[fips, 'median_width'], width/m2ft) / dfr.loc[fips, 'median_width'] * dfr.loc[fips, 'streetvalue_m_std']).sum() / 1e3 / valueRatio for width in widths ]
            nUnits = units.loc[fips, 'urbanized_units']
            stvalues_ph[fips] = [vv/nUnits*1e6 for vv in stvalues[fips]]
        stvalues['all'] = np.array([stvalues[ff] for ff in countylookup]).sum(axis=0)
        nUnits = units.loc[countylookup.keys(),'urbanized_units'].sum()
        stvalues_ph['all'] = [vv/nUnits*1e6 for vv in stvalues['all']]

        # sort plot order (and thus legend order) by the value at the max width
        valArray = np.array([stvalues_ph[fips][-1]*-1 for fips in countylookup])  # -1 to sort descending
        sortOrder= np.array(list(countylookup.keys()))[valArray.argsort()]

        fig, ax = plt.subplots(figsize=figsizes['full'])
        ax.plot(widths,  stvalues_ph['all'], lw=5, color='k', label='All counties')
        lstyles = ['solid','dotted','dashed','dashdot']*5
        colors = ['k']*7+['0.3']*7+['0.6']*7
        for ii, fips in enumerate(sortOrder):
            if self.bw:
                ax.plot(widths, stvalues_ph[fips], label=countylookup_formatted[fips], linestyle=lstyles[ii], color=colors[ii])
            else:
                ax.plot(widths, stvalues_ph[fips], label=countylookup_formatted[fips])
        ax.set_xlim(maxW,16) # reverse and truncate at 16 ft
        ax.set_ylim(0, 159.99)
        #ax.set_ylim(axlims)
        #ax2.set_ylim(0,axlims[1]*1e6/units.loc['all','urbanized_units'])
        ax.set_xlabel('Maximum residential street width (ft)', fontsize=12)
        #ax.set_ylabel('Total value (billion $)')
        ax.set_ylabel('Land value per housing unit\n(thousand $)', fontsize=12)
        ax.legend(bbox_to_anchor=(1.04,1.),loc='upper left', borderaxespad=0, fontsize=8, frameon=False)
        ax.text(-0.15, -0.2, 'x-axis is truncated at 16 feet, which is the approximate minimum required for access',ha='left', transform=ax.transAxes)
        fig.subplots_adjust(right=0.7, bottom=0.17)
        fig.savefig(paths['graphics']+'Fig6_valueofwidthreductions'+'_bw'*self.bw+'.pdf')

        # stats for text
        sc_value = dfr.loc['06085', 'streetvalue_m_std'].sum() * 1e3 / valueRatio/ units.loc['06085', 'urbanized_units']
        post2000_width = dfr[dfr.fbuy>=2000].loc['06085','median_width'].mean()
        hyp_value = (np.minimum(dfr.loc['06085', 'median_width'], post2000_width) / dfr.loc['06085', 'median_width'] * dfr.loc['06085', 'streetvalue_m_std']).sum() * 1e3 / valueRatio/ units.loc['06085', 'urbanized_units']
        print('Reducing Santa Clara streets to post-2000 width ({:.0f} ft) saves ${:.0f},000 per unit'.format(post2000_width*m2ft, sc_value-hyp_value))
        hyp_value = (np.minimum(dfr.loc['06085', 'median_width'], 16/m2ft) / dfr.loc['06085', 'median_width'] * dfr.loc['06085', 'streetvalue_m_std']).sum() * 1e3 / valueRatio/ units.loc['06085', 'urbanized_units']
        print('Reducing them to 16 ft saves ${:.0f},000 per unit'.format(sc_value-hyp_value))

    def landValueDistribution(self):
        """box plots of land value distribution within counties"""
        _ = self.getDf()
        dfr = self.dfr.copy()
        dfr['county'] = dfr.fips.map(countylookup_formatted)
        dfr['lv'] = dfr.landvalue_acre_std / 1e6

        fig, ax = plt.subplots(figsize=figsizes['full'])
        dfr.boxplot('lv', by='county', rot=90, grid=False, ax=ax)
        ax.set_xlabel('')
        ax.set_ylabel('Residential land value (million $ per acre)', fontsize=12)
        ax.set_title('')
        fig.suptitle('')
        fig.text(0.1,0.05,'Distribution is by street', fontsize=10,  style='italic')
        fig.tight_layout()
        fig.savefig(paths['graphics']+'FigA5_landvalue_distribution.pdf')

    def plotAll(self):
        self.widthHistograms()
        self.timeTrends()
        self.streetValue()
        self.landValueDistribution()

        # now do B&W versions for print journal
        self.bw=True
        print('Now doing black and white versions')
        self.timeTrends()
        self.streetValue()

def validation():
    """Get random sample of points in order to measure street width in Google Earth
    Validates the GIS-based estimates against hand-measurement in Google Earth"""
    outFn = paths['working']+'validation.csv'
    if os.path.exists(outFn):
        print('Validation file {} already exists. Delete it if you want to overwrite.'.format(outFn))
        return
        
    # Select random sample of points
    db = tools.dbConnection(user='adammb', db='streetwidths', host='localhost', schema='main', requirePassword=False)
    db.execute('SELECT setseed(0.5);')
    dfs = {}
    for fips in fipslookup.values():
        print('Getting sample for county {}'.format(fips))
        cmd = '''SELECT {fips} AS fips, p.osm_id, p.sid, p.frc, o.median_width, 
                    ST_X(ST_Transform(p.geom, 4326)) AS lon, ST_Y(ST_Transform(p.geom, 4326)) AS lat
                FROM streetpts_{fips} p LEFT JOIN osm_metrics_{fips} o USING (osm_id)
                WHERE ROUND(o.median_width)=ROUND(p.dist*2) AND o.clazz=41 AND urbanized
                ORDER BY RANDOM() LIMIT 10;'''.format(fips=fips)
        dfs[fips] = db.execfetchDf(cmd)

    df = pd.concat(dfs.values())
    df['latlong'] = df[['lat','lon']].astype(str).apply(','.join, axis=1)
    df['measured_width'] = ''
    df.to_csv(outFn, index=False)

def atlasChart(bw=False):
    """Bar chart of street widths based on Atlas of Urban Expansion data"""
    colors = ['0.1','0.3'] if bw else c2s
    colsToUse = ['Average Road Width (meters)','Unnamed: 11_level_0','City Name','Country','Region']
    df = pd.read_csv(paths['atlas']+'Blocks_and_Roads_Table_1.csv', header=[0,1])[colsToUse]
    df.columns = ['Width_pre1990','Width_1990-2015','city','country','region']
    df.Width_pre1990 = df.Width_pre1990 * m2ft
    df['Width_1990-2015'] = df['Width_1990-2015'] * m2ft
    df['us'] = df.country=='United States'
    df['color'] = df.us.map({True:colors[0],False:colors[1]})
    df.sort_values(by=['us','Width_1990-2015'], inplace=True)

    means = df.groupby('us')['Width_1990-2015'].mean()

    fig, ax = plt.subplots(figsize = figsizes['short'])
    #df[['Width_pre1990','Width_1990-2015']].plot.bar(rot=0)
    df['Width_1990-2015'].plot.bar(rot=0, ax = ax, color=df.color)
    ax.set_xticks([92,193])
    ax.tick_params(axis='x', length=0)
    ax.set_xticklabels(['Non-US cities','US'], ha='center', size=12)
    ax.set_ylabel('Mean ROW width (feet)', size=12)

    ax.hlines(means.loc[False], xmin=-0.5,xmax=200,lw=1, color=colors[1])
    ax.hlines(means.loc[True], xmin=-0.5,xmax=200,lw=1, color=colors[0])
    ax.text(10,25,'Non-US city mean', color=colors[1])
    ax.text(10,39.5,'US city mean', color=colors[0])

    fig.tight_layout()
    fig.savefig(paths['graphics']+'Fig1_atlas_widths'+'_bw'*bw+'.jpg', dpi=600, tight_layout=True)
    fig.savefig(paths['graphics']+'Fig1_atlas_widths'+'_bw'*bw+'.pdf', tight_layout=True)

def export_for_mapbox(fipslist = None):
    """Exports a GPKG for the Mapbox site"""
    if fipslist is None:
        fipslist = list(fipslookup.values())
    if isinstance(fipslist,str):
        fipslist = [fipslist]
    assert isinstance(fipslist, list)

    db = tools.dbConnection(user='adammb', db='streetwidths', host='localhost', schema='main', requirePassword=False)
    outPath = paths['website'] + 'Data/'
    for fips in fipslist:
        print('Exporting county {}'.format(fips))
        df = db.execfetchDf('SELECT osm_id, urbanized, ignore, fbuy, clazz, is_residential, length*{m2ft} AS length_ft, netlength*{m2ft} AS netlength_ft, width_ft, streetvalue_m_std, value_std_per_25ft FROM osm_metrics_{fips} WHERE ignore=False;'.format(m2ft=m2ft, fips=fips))
        df.width_ft = df.width_ft.round()
        df.to_csv(outPath+'osm_metrics_{}.tsv'.format(fips), sep='\t')
        cmd = '''ogr2ogr -f "GPKG" "{outPath}osm_{fips}.gpkg" "PG:dbname=streetwidths user=adammb" "rawdata.osm_{fips}"'''.format(outPath=outPath, fips=fips)
        assert os.system(cmd)==0    

def export_for_tableau():
    """Gets all the tsvs for Mapbox, merges them, and adds a county name field
    Note that the value histogram is exported in analysis()"""
    outPath = paths['website'] + 'Data/'
    fipslist = list(fipslookup.values())
    dfs = []
    for fips in fipslist:
        df = pd.read_csv(outPath+'osm_metrics_{}.tsv'.format(fips), sep='\t', index_col=None)
        df['county'] = countylookup_formatted[fips].replace(',',' County,')
        df['countyfips'] = fips
        dfs.append(df)
    df = pd.concat(dfs)
    df = df[df.urbanized]
    df.fbuy.fillna(-1, inplace=True)
    #df.drop(columns=['Unnamed: 0','osm_id','ignore', 'urbanized','clazz',inplace=True)
    # group by county, fbuy and is_residential to speed rendering
    aggdf = df.groupby(['county','countyfips','fbuy','is_residential','width_ft'])[['length_ft','streetvalue_m_std']].sum().reset_index()
    aggdf.to_csv(outPath+'metrics_aggregated.tsv', sep='\t', index=False)



if __name__ == '__main__':
   dataLoader().loadAll()
   createStreetPolygons().doAllCounties()
   analysis().plotAll()
   validation()
   atlasChart()
   atlasChart(bw=True)

