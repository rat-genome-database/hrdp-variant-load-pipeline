package edu.mcw.rgd.hrdpVariantLoad;

import edu.mcw.rgd.dao.DataSourceFactory;
import edu.mcw.rgd.datamodel.RgdId;
import edu.mcw.rgd.datamodel.Sample;
import edu.mcw.rgd.datamodel.variants.VariantMapData;
import edu.mcw.rgd.datamodel.variants.VariantSampleDetail;
import edu.mcw.rgd.process.Utils;
import edu.mcw.rgd.util.Zygosity;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.zip.GZIPInputStream;

public class HrdpVariants {

    private String version;
    int sampleIdStart;
    private int mapKey;
    protected Logger logger = LogManager.getLogger("status");
    private Zygosity zygosity = new Zygosity();
    private String inputDir;
    private String inputFile;
    private Map<String, Integer> colNameToSampleId;
    private Map<String, Integer> colNameToSampleId7;
    private final DAO dao = new DAO();
    private SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    private int zeroDepthCnt = 0;
    private int totalSamples = 0;
    public void main() throws Exception {

        logger.info(getVersion());
        logger.info("   "+dao.getConnection());

        long pipeStart = System.currentTimeMillis();
        logger.info("Pipeline started at "+sdt.format(new Date(pipeStart))+"\n");

        geneCacheMap = new HashMap<>();
        // loops through files
        File file = new File(inputFile);
        parse(file);
//        ArrayList<File> files = new ArrayList<>();
//        dao.listFilesInFolder(folder, files);
//
//        for (File file : files) {
//            parse(file);
//        }

        logger.info("Total pipeline runtime -- elapsed time: "+
                Utils.formatElapsedTime(pipeStart,System.currentTimeMillis()));
    }

    void parse(File file) throws Exception{
        zeroDepthCnt = 0;
        totalSamples = 0;
        // create samples and insert them
        // parse through files and insert variants
//        String name = getStrainName(file.getName());
//        System.out.println(name);
//        Integer strainRgdId = getStrainRgdId(name);
//
//        Sample sample = dao.getSampleByAnalysisNameAndMapKey(name,mapKey);
//
//        if (sample == null) {
//            logger.info("\t\tCreating new Sample for " + name);
//            sample = new Sample();
//            sample.setId(sampleIdStart);
//            sampleIdStart++;
//            sample.setAnalysisName(name);
//
//            if (strainRgdId != 0) {
//                sample.setStrainRgdId(strainRgdId);
//            }
//            sample.setGender("U");
//            sample.setDescription("Dr. Mindy Dwinell - Hybrid rat diversity program");
//            sample.setPatientId(380); // default is rat 8 patient id, 380
//            sample.setMapKey(mapKey);
//            sample.setGrantNumber("R24OD022617");
//            dao.insertSample(sample);
//        }

        // parse file and insert/add samples to variants
        logger.info("\tBegin parsing file... " + file.getName());
        BufferedReader br = dao.openFile(file.getAbsolutePath());
        String lineData;
        List<VariantMapData> variants = new ArrayList<>();
        List<VariantMapData> tobeUpdated = new ArrayList<>();
        List<VariantSampleDetail> samples = new ArrayList<>();
        HashMap<Integer, Sample> strainSamples = new HashMap<>();
        int totalVars = 0;
        while ((lineData = br.readLine()) != null){
            if (lineData.startsWith("##"))
                continue;
            if (lineData.startsWith("#CHROM")){
                // get sample from column names, and assign to map
                String[] splitCols = lineData.split("\t");
                for (int i = 9; i < splitCols.length; i++){
                    Integer sampleId = colNameToSampleId.get(splitCols[i]);
                    if (sampleId != null) {
                        Sample s = dao.getSampleBySampleId(sampleId);
                        strainSamples.put(i,s);
                    }
                    else
                        System.out.println(splitCols[i] + " not found!");
                }
                continue;
            }

            totalVars += parseLineData(lineData, strainSamples, tobeUpdated);


        } // end of file read
        if (zeroDepthCnt!=0){
            logger.info("\t\t\tVariants with 0 depth being ignored: "+zeroDepthCnt);
        }
        if (!tobeUpdated.isEmpty()){
            logger.info("\t\tVariants end pos being updated: "+tobeUpdated.size());
//            dao.updateVariantEndPosBatch(tobeUpdated);
        }
        if (totalVars != 0){
            logger.info("\t\tVariants being entered: " + totalVars);
//            dao.insertVariantRgdIds(variants);
//            dao.insertVariants(variants);
//            dao.insertVariantMapData(variants);
        }
        if (totalSamples != 0){
//            totalSamples += samples.size();
            logger.info("\t\tTotal samples being created: " + totalSamples);
//            dao.insertVariantSample(samples);
        }
        logger.info("\tEnd parsing file... " + file.getName());
//        System.out.println(variants.size());
    }

    String getStrainName(String fileName) throws Exception{
        String strain = fileName;
        if (strain.contains("_PASS"))
        {
            strain = strain.replace("_PASS","");
        }
        int lastUnderScore = strain.lastIndexOf('_');
        strain = strain.substring(0,lastUnderScore);
        lastUnderScore = strain.lastIndexOf('_');
        strain = strain.substring(0,lastUnderScore)+")";
        int cnt =0;
        for (int i = 0; i < strain.length(); i++){
            if (strain.charAt(i) == '_'){
                cnt++;
            }
        }
        if (cnt>2) {
            strain = strain.replaceFirst("_", "-");
        }
        strain = strain.replaceFirst("_","/");
        strain = strain.replace("_"," (");

        return strain;
    }

    Integer parseLineData(String lineData, HashMap<Integer,Sample> sampleMap, List<VariantMapData> tobeUpdated) throws Exception{
        // CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ACI_EurMcwi_2019NG
        VariantMapData v = new VariantMapData();
        ArrayList<VariantMapData> vars = new ArrayList<>();
        ArrayList<VariantMapData> existing = new ArrayList<>();
        List<VariantSampleDetail> samples = new ArrayList<>();
        boolean needCopyVar = false;
        boolean needCopyRef = false;
        Integer totalDepth = null;
        String[] data = lineData.split("\t");
        String[] depths = new String[0];
        // first check if ref or alt has a ','
        // if yes, make a copy of
        for (int i = 0; i < 9; i++){
            switch (i){
                case 0: // chrom
                    // chrM is for MT
                    if (data[i].contains("unplaced") || data[i].contains("unloc") || data[i].contains("contig") || data[i].contains("scaffold")){
                        return 0;
                    }
                    v.setChromosome(data[i].replace("chr",""));
                    if (v.getChromosome().equalsIgnoreCase("M"))
                        v.setChromosome("MT");
                    break;
                case 1: // pos
                    int start = Integer.parseInt(data[i]);
                    v.setStartPos(start);
                    break;
                case 2: // id
                    if (!data[i].equals(".")) {
                        v.setRsId(data[i]);
                    }
                    break;
                case 3: // ref
                    if (data[i].contains(",")) {
                        needCopyRef = true;
                        v.setReferenceNucleotide(data[i]); // change in alt for padding base if it exists
                    }
                    else {
                        v.setReferenceNucleotide(data[i]); // change in alt for padding base if it exists
                    }
                    break;
                case 4: // alt
                    // get end position, insertion - start+1, deletion/SNP - length of alt/var
                    // also set padding base if possible
                    // make sure to copy variant if there is a ','
                    // '*' represents deletion
                    if (data[i].contains(",")) {
                        needCopyVar = true;
                        v.setVariantNucleotide(data[i]);
                    }
                    else {

                        if (data[i].equals("*")) {
                            v.setVariantNucleotide(null);
                            v.setEndPos(v.getStartPos() + v.getReferenceNucleotide().length());
                            v.setVariantType("deletion");
                        }
                        else {
                            String var = data[i];
                            v.setVariantNucleotide(var);
                            if (v.getReferenceNucleotide().length() > var.length() && var.length() == 1) {
                                // deletion
                                v.setPaddingBase(var);
                                v.setStartPos(v.getStartPos()+1);
                                v.setVariantNucleotide(null);
                                String ref = v.getReferenceNucleotide().substring(1);
                                v.setReferenceNucleotide(ref);
                                v.setEndPos(v.getStartPos() + v.getReferenceNucleotide().length());
                                v.setVariantType("deletion");
                            } else if (v.getReferenceNucleotide().length() > var.length() && v.getReferenceNucleotide().startsWith(var)){
                                v.setPaddingBase(var);
                                v.setStartPos(v.getStartPos()+var.length());
                                v.setVariantNucleotide(null);
                                String ref = v.getReferenceNucleotide().replaceFirst(var, "");
                                v.setReferenceNucleotide(ref);
                                v.setEndPos(v.getStartPos()+ref.length());
                                v.setVariantType("deletion");
                            } else if (var.length() > v.getReferenceNucleotide().length() && v.getReferenceNucleotide().length() == 1) {
                                // insertion
                                v.setPaddingBase(v.getReferenceNucleotide());
                                v.setStartPos(v.getStartPos()+1);
                                v.setEndPos(v.getStartPos()+1);
                                v.setReferenceNucleotide(null);
                                var = var.substring(1);
                                v.setVariantNucleotide(var);
                                v.setVariantType("insertion");
                            } else if (var.length() > v.getReferenceNucleotide().length() && var.startsWith(v.getReferenceNucleotide())){
                                v.setPaddingBase(v.getReferenceNucleotide());
                                v.setStartPos(v.getStartPos()+v.getReferenceNucleotide().length());
                                v.setEndPos(v.getStartPos()+1);
                                v.setReferenceNucleotide(null);
                                var = var.replaceFirst(v.getPaddingBase(),"");
                                v.setVariantNucleotide(var);
                                v.setVariantType("insertion");
                            } else {
                                v.setVariantNucleotide(data[i]);
                                if (v.getReferenceNucleotide().length() == v.getVariantNucleotide().length()) {
                                    v.setEndPos(v.getStartPos() + 1);
                                    if (v.getReferenceNucleotide().length() > 1) {
                                        v.setVariantType("mnv");
                                        v.setEndPos(v.getStartPos()+v.getReferenceNucleotide().length());
                                    }
                                    else {
                                        v.setVariantType("snv");
                                        v.setEndPos(v.getStartPos()+1);
                                    }
                                } else if (v.getReferenceNucleotide().length() > v.getVariantNucleotide().length()) {
                                    v.setEndPos(v.getStartPos() + v.getReferenceNucleotide().length());
                                    v.setVariantType("delins");
                                } else {
                                    v.setEndPos(v.getStartPos() + 1);
                                    v.setVariantType("delins");
                                }
                            }
                        }
                    }
                    break;
                case 5: // qual
                case 6: // filter
                case 7: // info
                    // In the INFO field (long line of annotations added in the analysis)  you will find AF – allele frequency 0.5 it is HET
                    // AC=1; AF=0.500 ;AN=2;ANN=T
                case 8: //format
                    // In the FORMAT field you will find DP – reads depth (total number of reads per this position)
                    break;
                case 9: // actual format data
                    // 0/1 :ref 32,allele 9: total depth 41 :99:130,0,872

                    String[] formatData = data[i].split(":");
                    depths = formatData[1].split(","); // first is ref, following are alleles
                    totalDepth = Integer.parseInt(formatData[2]);
                    if (totalDepth==0)
                    {
                        zeroDepthCnt++;
//                        logger.info(lineData);
                        return 0;
                    }
                    break;
            }
        }
        if (isGenic(v))
            v.setGenicStatus( "GENIC");
        else
            v.setGenicStatus("INTERGENIC" );
        v.setMapKey(mapKey);
        v.setSpeciesTypeKey(3);
        List<VariantMapData> dbVars;
        if (Utils.isStringEmpty(v.getRsId()))
             dbVars = dao.getVariant(v);
        else
            dbVars = dao.getVariantByRsId(v);
        List<VariantMapData> newVars = new ArrayList<>();
        if (needCopyRef || needCopyVar){
            String[] refNucs = {};
            String[] varNucs = {};
            int cnt;
            if (needCopyRef) {
                refNucs = v.getReferenceNucleotide().split(",");
                cnt = refNucs.length;
            }
            else {
                varNucs = v.getVariantNucleotide().split(",");
                cnt = varNucs.length;
            }


//            List<VariantMapData> variantCopies = new ArrayList<>();
            for (int i = 0; i < cnt; i++){
                VariantMapData copy = new VariantMapData();
                copy.setChromosome(v.getChromosome());
                copy.setRsId(v.getRsId());

                copy.setStartPos(v.getStartPos());
                String var = null;
                if (needCopyRef){
                    copy.setReferenceNucleotide(refNucs[i]);
                    var = v.getVariantNucleotide();
                }
                else{
                    copy.setReferenceNucleotide(v.getReferenceNucleotide());
                    var = varNucs[i];
                }
                if (var.equals("*")){
                    copy.setVariantNucleotide(null);
                    copy.setEndPos(copy.getStartPos() + copy.getReferenceNucleotide().length());
                    copy.setVariantType("deletion");
                }
                else {

                    if (copy.getReferenceNucleotide().length() > var.length() && var.length() == 1) {
                        // deletion
                        copy.setPaddingBase(var);
                        copy.setVariantNucleotide(null);
                        String ref = copy.getReferenceNucleotide().substring(1);
                        copy.setReferenceNucleotide(ref);
                        copy.setEndPos(copy.getStartPos() + copy.getReferenceNucleotide().length());
                        copy.setVariantType("deletion");
                    } else if (v.getReferenceNucleotide().length() > var.length() && v.getReferenceNucleotide().startsWith(var)){
                        copy.setPaddingBase(var);
                        copy.setStartPos(v.getStartPos()+var.length());
                        copy.setVariantNucleotide(null);
                        String ref = v.getReferenceNucleotide().replaceFirst(var, "");
                        copy.setReferenceNucleotide(ref);
                        copy.setEndPos(v.getStartPos()+ref.length());
                        copy.setVariantType("deletion");
                    } else if (var.length() > copy.getReferenceNucleotide().length() && copy.getReferenceNucleotide().length() == 1) {
                        // insertion
                        copy.setPaddingBase(v.getReferenceNucleotide());
                        copy.setEndPos(v.getStartPos()+1);
                        copy.setReferenceNucleotide(null);
                        var = var.substring(1);
                        copy.setVariantNucleotide(var);
                        copy.setVariantType("insertion");
                    } else if (var.length() > v.getReferenceNucleotide().length() && var.startsWith(v.getReferenceNucleotide())){
                        copy.setPaddingBase(v.getReferenceNucleotide());
                        copy.setStartPos(v.getStartPos()+v.getReferenceNucleotide().length());
                        copy.setEndPos(v.getStartPos()+1);
                        copy.setReferenceNucleotide(null);
                        var = var.replaceFirst(copy.getPaddingBase(),"");
                        copy.setVariantNucleotide(var);
                        copy.setVariantType("insertion");
                    } else {
                        copy.setVariantNucleotide(varNucs[i]);
                        if (copy.getReferenceNucleotide().length() == copy.getVariantNucleotide().length()) {
                            copy.setEndPos(copy.getStartPos() + 1);
                            if (copy.getReferenceNucleotide().length() > 1) {
                                copy.setEndPos(copy.getStartPos()+copy.getReferenceNucleotide().length());
                                copy.setVariantType("mnv");
                            }
                            else {
                                copy.setEndPos(v.getStartPos()+1);
                                copy.setVariantType("snp");
                            }
                        } else if (copy.getReferenceNucleotide().length() > copy.getVariantNucleotide().length()) {
                            copy.setEndPos(copy.getStartPos() + copy.getReferenceNucleotide().length());
                            copy.setVariantType("delins");
                        } else {
                            copy.setEndPos(copy.getStartPos() + 1);
                            copy.setVariantType("delins");
                        }
                    }
                }
                copy.setMapKey(mapKey);
                copy.setGenicStatus(v.getGenicStatus());
                copy.setSpeciesTypeKey(3);
//                variantCopies.add(copy);
                boolean exist = false;
                for (VariantMapData dbVar : dbVars) {
                    if (Utils.stringsAreEqual(copy.getReferenceNucleotide(),dbVar.getReferenceNucleotide() ) && Utils.stringsAreEqual(copy.getVariantNucleotide(), dbVar.getVariantNucleotide())
                    && copy.getStartPos()==dbVar.getStartPos()){
                        exist = true;
                        existing.add(dbVar);
                        if (dbVar.getEndPos() != copy.getEndPos() && copy.getEndPos()!=0) {
                            dbVar.setEndPos(copy.getEndPos());
                            tobeUpdated.add(dbVar);
                        }
                        break;
                    }
                }// end check with database vars
                if (!exist) {
                    RgdId r = dao.createRgdId(RgdId.OBJECT_KEY_VARIANTS, "ACTIVE", "created by HRDP Load Pipeline", mapKey);
                    copy.setId(r.getRgdId());
//                    v.setId(0);
                    newVars.add(copy);

                }

            }

            vars.addAll(newVars);
        }
        else {
                boolean exist = false;
                for (VariantMapData dbVar : dbVars) {
                    if (Utils.stringsAreEqual(v.getReferenceNucleotide(), dbVar.getReferenceNucleotide()) && Utils.stringsAreEqual(v.getVariantNucleotide(), dbVar.getVariantNucleotide())
                            && v.getStartPos() == dbVar.getStartPos()) {
                        exist = true;
                        existing.add(dbVar);
//
//                        }
                        if (dbVar.getEndPos() != v.getEndPos() && v.getEndPos() != 0) {
                            dbVar.setEndPos(v.getEndPos());
                            tobeUpdated.add(dbVar);
                        }
                        break;
                    }
                }
                if (!exist) {

                    RgdId r = dao.createRgdId(RgdId.OBJECT_KEY_VARIANTS, "ACTIVE", "created by HRDP Load Pipeline", mapKey);
                    v.setId(r.getRgdId());
//                    v.setId(0);
                    vars.add(v);

                }
//            vars.add(v);

        }
        ArrayList<VariantMapData> varsForSamples = new ArrayList<>();
        varsForSamples.addAll(vars);
        varsForSamples.addAll(existing);
        for (int i = 9; i < data.length; i++) {
            String[] formatData = data[i].split(":");
            if (Utils.stringsAreEqual(formatData[0],"0/0") || Utils.stringsAreEqual(formatData[0],"./."))
                continue;
            depths = formatData[1].split(","); // first is ref, following are alleles
            try {
                totalDepth = Integer.parseInt(formatData[2]);
            }
            catch (Exception ignore){
                System.out.println("");
            } // depth is "."
            Sample s = sampleMap.get(i);
            for (int j = 0; j < varsForSamples.size(); j++) {
                VariantMapData dbVar = varsForSamples.get(j);
                int varFreq = Integer.parseInt(depths[j+1]);
                int sampleCnt = dao.getVariantSampleDetailCount((int) dbVar.getId(), s.getId());
                if (sampleCnt == 0 && varFreq != 0) {
                    VariantSampleDetail newSample = new VariantSampleDetail();
                    newSample.setId(dbVar.getId());
                    newSample.setSampleId(s.getId());
                    newSample.setDepth(totalDepth);
                    newSample.setVariantFrequency(varFreq);
                    zygosity.computeZygosityStatus(newSample.getVariantFrequency(), newSample.getDepth(), s.getGender(), dbVar, newSample);

                    int zygPercentRead = varFreq / newSample.getDepth();
                    newSample.setZygosityPercentRead(zygPercentRead);
                    samples.add(newSample);

                }
            }
        }
        if (!vars.isEmpty()){
            dao.insertVariantRgdIds(vars);
            dao.insertVariants(vars);
            dao.insertVariantMapData(vars);
        }
        if (!samples.isEmpty()){
            totalSamples += samples.size();
            dao.insertVariantSample(samples);
        }
        return vars.size();
    }

    Integer getStrainRgdId(String sampleName) throws Exception {
        int strainStop = sampleName.indexOf('(')-1;
        String strainName = sampleName.substring(0,strainStop);
        return dao.getStrainRgdIdByTaglessStrainSymbol(strainName);
    }

    boolean isGenic(VariantMapData v) throws Exception {

        GeneCache geneCache = geneCacheMap.get(v.getChromosome());
        if( geneCache==null ) {
            geneCache = new GeneCache();
            geneCacheMap.put(v.getChromosome(), geneCache);
            geneCache.loadCache(mapKey, v.getChromosome(), DataSourceFactory.getInstance().getDataSource());
        }
        List<Integer> geneRgdIds = geneCache.getGeneRgdIds((int)v.getStartPos(), (int)v.getEndPos());
        logger.debug("Variant position: "+v.getChromosome()+":"+v.getStartPos());
        for (int id : geneRgdIds){
            logger.debug("\tGene rgdId: "+id);
        }
        return !geneRgdIds.isEmpty();
    }
    Map<String, GeneCache> geneCacheMap;

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }

    public void setSampleIdStart(int sampleIdStart) {
        this.sampleIdStart = sampleIdStart;
    }

    public int getSampleIdStart() {
        return sampleIdStart;
    }

    public void setInputDir(String inputDir) {
        this.inputDir = inputDir;
    }

    public String getInputDir(){
        return inputDir;
    }

    public void setMapKey(int mapKey) {
        this.mapKey = mapKey;
    }

    public int getMapKey() {
        return mapKey;
    }

    public void setInputFile(String inputFile) {
        this.inputFile = inputFile;
    }

    public String getInputFile() {
        return inputFile;
    }

    public Map<String, Integer> getColNameToSampleId() {
        return colNameToSampleId;
    }

    public void setColNameToSampleId(Map<String, Integer> colNameToSampleId) {
        this.colNameToSampleId = colNameToSampleId;
    }

    public Map<String, Integer> getColNameToSampleId7() {
        return colNameToSampleId7;
    }

    public void setColNameToSampleId7(Map<String, Integer> colNameToSampleId7) {
        this.colNameToSampleId7 = colNameToSampleId7;
    }
}
