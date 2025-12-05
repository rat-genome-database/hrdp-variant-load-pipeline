package edu.mcw.rgd.hrdpVariantLoad;

import edu.mcw.rgd.dao.DataSourceFactory;
import edu.mcw.rgd.dao.impl.*;
import edu.mcw.rgd.dao.impl.variants.VariantDAO;
import edu.mcw.rgd.dao.spring.variants.VariantMapQuery;
import edu.mcw.rgd.dao.spring.variants.VariantSampleQuery;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.datamodel.variants.VariantMapData;
import edu.mcw.rgd.datamodel.variants.VariantSampleDetail;
import org.springframework.jdbc.core.SqlParameter;
import org.springframework.jdbc.object.BatchSqlUpdate;

import javax.sql.DataSource;
import java.io.*;
import java.sql.Types;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.zip.GZIPInputStream;

public class DAO {
    private SampleDAO sampleDAO = new SampleDAO();
    private StrainDAO sdao = new StrainDAO();
    private VariantDAO vdao = new VariantDAO();
    private RGDManagementDAO managementDAO = new RGDManagementDAO();
    private MapDAO mdao = new MapDAO();
    private GeneDAO gdao = new GeneDAO();
    public String getConnection() {
        return vdao.getConnectionInfo();
    }

    public DataSource getVariantDataSource() throws Exception{
        return DataSourceFactory.getInstance().getCarpeNovoDataSource();
    }
    public Sample getSampleByAnalysisNameAndMapKey(String name, int mapKey) throws Exception{
        sampleDAO.setDataSource(getVariantDataSource());
        return sampleDAO.getSampleByAnalysisNameAndMapKey(name,mapKey);
    }

    public Sample getSampleBySampleId(int sampleId) throws Exception {
        sampleDAO.setDataSource(getVariantDataSource());
        return sampleDAO.getSampleBySampleId(sampleId);
    }

    public Integer getStrainRgdIdByTaglessStrainSymbol(String strainName) throws Exception {
        return sdao.getStrainRgdIdByTaglessStrainSymbol(strainName);
    }

    public void insertSample(Sample sample) throws Exception {
        vdao.insertSample(sample);
    }

    public RgdId createRgdId(int objectKey, String objectStatus, String notes, int mapKey) throws Exception{
        int speciesKey = SpeciesType.getSpeciesTypeKeyForMap(mapKey);
        return managementDAO.createRgdId(objectKey, objectStatus, notes, speciesKey);
    }

    public VariantSampleDetail getVariantSampleDetail(int rgdId, int sampleId) throws Exception{
        return vdao.getVariantSampleDetailByRGDIdSampleId(rgdId,sampleId);
    }

    public int getVariantSampleDetailCount(int rgdId, int sampleId) throws Exception{
        return vdao.getVariantSampleDetailCount(rgdId, sampleId);
    }

    public int insertVariantSample(List<VariantSampleDetail> sampleData) throws Exception {
        BatchSqlUpdate bsu= new BatchSqlUpdate(this.getVariantDataSource(),
                "INSERT INTO variant_sample_detail (\n" +
                        " RGD_ID,SOURCE,SAMPLE_ID,TOTAL_DEPTH,VAR_FREQ,ZYGOSITY_STATUS,ZYGOSITY_PERCENT_READ," +
                        "ZYGOSITY_POSS_ERROR,ZYGOSITY_REF_ALLELE,ZYGOSITY_NUM_ALLELE,ZYGOSITY_IN_PSEUDO,QUALITY_SCORE)\n" +
                        "VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                new int[]{Types.INTEGER,Types.VARCHAR,Types.INTEGER, Types.INTEGER, Types.INTEGER,Types.VARCHAR, Types.INTEGER,
                        Types.VARCHAR,Types.VARCHAR, Types.INTEGER,Types.VARCHAR, Types.INTEGER}, 10000);
        bsu.compile();
        for(VariantSampleDetail v: sampleData ) {
            bsu.update(v.getId(), v.getSource(), v.getSampleId(),v.getDepth(),v.getVariantFrequency(),v.getZygosityStatus(),v.getZygosityPercentRead(),
                    v.getZygosityPossibleError(),v.getZygosityRefAllele(),v.getZygosityNumberAllele(),v.getZygosityInPseudo(),v.getQualityScore());
        }
        bsu.flush();
        // compute nr of rows affected
        int totalRowsAffected = 0;
        for( int rowsAffected: bsu.getRowsAffected() ) {
            totalRowsAffected += rowsAffected;
        }
        return totalRowsAffected;
    }

    public void insertVariantRgdIds(Collection<VariantMapData> mapData) throws Exception {
        vdao.insertVariantRgdIds(mapData);
    }
    public void insertVariants(List<VariantMapData> mapsData)  throws Exception{
        BatchSqlUpdate sql1 = new BatchSqlUpdate(this.getVariantDataSource(),
                "INSERT INTO variant (\n" +
                        " RGD_ID,REF_NUC, VARIANT_TYPE, VAR_NUC, RS_ID, CLINVAR_ID, SPECIES_TYPE_KEY)\n" +
                        "VALUES (?,?,?,?,?,?,?)",
                new int[]{Types.INTEGER,Types.VARCHAR,Types.VARCHAR, Types.VARCHAR, Types.VARCHAR, Types.VARCHAR,Types.INTEGER}, 10000);
        sql1.compile();
        for( VariantMapData v: mapsData) {
            long id = v.getId();
            sql1.update(id, v.getReferenceNucleotide(), v.getVariantType(), v.getVariantNucleotide(), v.getRsId(), v.getClinvarId(), v.getSpeciesTypeKey());

        }
        sql1.flush();
    }
    public void insertVariantMapData(List<VariantMapData> mapsData)  throws Exception{
        BatchSqlUpdate sql2 = new BatchSqlUpdate(this.getVariantDataSource(),
                "INSERT INTO variant_map_data (\n" +
                        " RGD_ID,CHROMOSOME,START_POS,END_POS,PADDING_BASE,GENIC_STATUS,MAP_KEY)\n" +
                        "VALUES (?,?,?,?,?,?,?)",
                new int[]{Types.INTEGER,Types.VARCHAR, Types.INTEGER, Types.INTEGER, Types.VARCHAR,Types.VARCHAR, Types.INTEGER}, 10000);
        sql2.compile();
        for( VariantMapData v: mapsData) {
            long id = v.getId();
            sql2.update(id, v.getChromosome(), v.getStartPos(), v.getEndPos(), v.getPaddingBase(), v.getGenicStatus(), v.getMapKey());
        }
        sql2.flush();
    }

    public List<VariantMapData> getVariant(VariantMapData v)throws Exception{
        String sql = "SELECT * FROM variant v inner join variant_map_data vmd on v.rgd_id=vmd.rgd_id where vmd.map_key=? and vmd.chromosome=? and vmd.start_pos=?";
        VariantMapQuery q = new VariantMapQuery(getVariantDataSource(), sql);
        q.declareParameter(new SqlParameter(Types.INTEGER));
        q.declareParameter(new SqlParameter(Types.VARCHAR));
        q.declareParameter(new SqlParameter(Types.INTEGER));
        return q.execute(v.getMapKey(), v.getChromosome(), v.getStartPos());
    }

    public List<VariantMapData> getVariantByRsId(VariantMapData v) throws Exception {
        String sql = "SELECT * FROM variant v inner join variant_map_data vmd on v.rgd_id=vmd.rgd_id where vmd.map_key=? v.rs_id=?";
        VariantMapQuery q = new VariantMapQuery(getVariantDataSource(), sql);
        q.declareParameter(new SqlParameter(Types.INTEGER));
        q.declareParameter(new SqlParameter(Types.VARCHAR));
        return q.execute(v.getMapKey(), v.getRsId());
    }

    public List<VariantMapData> getVariantsWithGeneLocation(int mapKey, String chr, int start, int stop) throws Exception {
        return vdao.getVariantsWithGeneLocation(mapKey,chr, start, stop);
    }

    public void updateVariantEndPosBatch(Collection<VariantMapData> toBeUpdated) throws Exception{
        BatchSqlUpdate su = new BatchSqlUpdate(this.getVariantDataSource(),
                "update variant_map_data set END_POS=? where RGD_ID=?",
                new int[]{Types.INTEGER, Types.INTEGER}, 5000);
        su.compile();
        for (VariantMapData v : toBeUpdated){
            su.update(v.getEndPos(),v.getId());
        }
        su.flush();
    }

    public void updateVariantMapDataGenicStatus(List<VariantMapData> mapsData) throws Exception {
        BatchSqlUpdate sql2 = new BatchSqlUpdate(this.getVariantDataSource(),
                "update variant_map_data set GENIC_STATUS=? where RGD_ID=?",
                new int[]{Types.VARCHAR,Types.INTEGER});
        sql2.compile();
        for( VariantMapData v: mapsData) {
            long id = v.getId();
            sql2.update(v.getGenicStatus(),id);
        }
        sql2.flush();
    }

    public List<MapData> getMapDataWithinRange(int start, int stop, String chr, int mapKey, int range) throws Exception{
        return mdao.getMapDataWithinRange(start,stop,chr,mapKey,range);
    }

    public Gene getGene(int rgdId) throws Exception{
        return gdao.getGene(rgdId);
    }

    public void listFilesInFolder(File folder, ArrayList<File> vcfFiles) throws Exception {
        for (File file : Objects.requireNonNull(folder.listFiles())) {
            if (file.isDirectory()) {
                listFilesInFolder(file,vcfFiles);
            } else {
                if (file.getName().endsWith(".vcf.gz")) {
//                    System.out.println(file.getName());
                    vcfFiles.add(file);
                }
            }
        }
    }

    public BufferedReader openFile(String fileName) throws IOException {

        String encoding = "UTF-8"; // default encoding

        InputStream is;
        if( fileName.endsWith(".gz") ) {
            is = new GZIPInputStream(new FileInputStream(fileName));
        } else {
            is = new FileInputStream(fileName);
        }

        BufferedReader reader = new BufferedReader(new InputStreamReader(is, encoding));
        return reader;
    }
}
