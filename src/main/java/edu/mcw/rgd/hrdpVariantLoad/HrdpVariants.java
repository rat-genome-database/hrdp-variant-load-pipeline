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
    protected Logger logger = LogManager.getLogger("status");
    private Zygosity zygosity = new Zygosity();
    private String inputDir;
    private int mapKey = 0;
    private final DAO dao = new DAO();
    private SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    private int zeroDepthCnt = 0;
    public void main(String[] args) throws Exception {

        logger.info(getVersion());
        logger.info("   "+dao.getConnection());

        long pipeStart = System.currentTimeMillis();
        logger.info("Pipeline started at "+sdt.format(new Date(pipeStart))+"\n");

        for (int i = 0; i < args.length; i++){
            switch (args[i]){
                case "--mapKey":
                    mapKey = Integer.parseInt(args[++i]);
                    break;
            }
        }

        geneCacheMap = new HashMap<>();
        // loops through files
        File folder = new File(inputDir);
        ArrayList<File> files = new ArrayList<>();
        listFilesInFolder(folder, files);

        for (File file : files) {
            parse(file);
        }

        logger.info("Total pipeline runtime -- elapsed time: "+
                Utils.formatElapsedTime(pipeStart,System.currentTimeMillis()));
    }

    void parse(File file) throws Exception{
        zeroDepthCnt = 0;
        // create samples and insert them
        // parse through files and insert variants
        String name = getStrainName(file.getName());
//        System.out.println(name);
        Integer strainRgdId = getStrainRgdId(name);

        Sample sample = dao.getSampleByAnalysisNameAndMapKey(name,mapKey);

        if (sample == null) {
            logger.info("\t\tCreating new Sample for " + name);
            sample = new Sample();
            sample.setId(sampleIdStart);
            sampleIdStart++;
            sample.setAnalysisName(name);

            if (strainRgdId != 0) {
                sample.setStrainRgdId(strainRgdId);
            }
            sample.setGender("U");
            sample.setDescription("Dr. Mindy Dwinell - Hybrid rat diversity program");
            sample.setPatientId(372); // default is rat 7.2 patient id, 372
            sample.setMapKey(mapKey);
            sample.setGrantNumber("R24OD022617");
            dao.insertSample(sample);
        }

        // parse file and insert/add samples to variants
        logger.info("\tBegin parsing file... " + file.getName());
        BufferedReader br = openFile(file.getAbsolutePath());
        String lineData;
        List<VariantMapData> variants = new ArrayList<>();
        List<VariantMapData> tobeUpdated = new ArrayList<>();
        List<VariantSampleDetail> samples = new ArrayList<>();
        while ((lineData = br.readLine()) != null){
            if (lineData.startsWith("#"))
                continue;


            List<VariantMapData> vars = parseLineData(lineData, samples, sample, tobeUpdated);
            if (!vars.isEmpty())
                variants.addAll(vars);

        } // end of file read
        if (zeroDepthCnt!=0){
            logger.info("\t\t\tVariants with 0 depth being ignored: "+zeroDepthCnt);
        }
        if (!tobeUpdated.isEmpty()){
            logger.info("\t\tVariants end pos being updated: "+tobeUpdated.size());
            dao.updateVariantEndPosBatch(tobeUpdated);
        }
        if (!variants.isEmpty()){
            logger.info("\t\tVariants being entered: " + variants.size());
            dao.insertVariants(variants);
            dao.insertVariantMapData(variants);
        }
        if (!samples.isEmpty()){
            logger.info("\t\tNew samples being created: " + samples.size());
            dao.insertVariantSample(samples);
        }
        logger.info("\tEnd parsing file... " + file.getName());
//        System.out.println(variants.size());
    }

    void listFilesInFolder(File folder, ArrayList<File> vcfFiles) throws Exception {
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

    List<VariantMapData> parseLineData(String lineData, List<VariantSampleDetail> samples, Sample s, List<VariantMapData> tobeUpdated) throws Exception{
        // CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ACI_EurMcwi_2019NG
        VariantMapData v = new VariantMapData();
        List<VariantMapData> vars = new ArrayList<>();
        boolean needCopyVar = false;
        boolean needCopyRef = false;
        Integer totalDepth = null;
        String[] data = lineData.split("\t");
        String[] depths = new String[0];
        // first check if ref or alt has a ','
        // if yes, make a copy of
        for (int i = 0; i < data.length; i++){
            switch (i){
                case 0: // chrom
                    // chrM is for MT
                    if (data[i].contains("unplaced") || data[i].contains("unloc") || data[i].contains("contig") || data[i].contains("scaffold")){
                        return new ArrayList<>(0);
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
                    if (!data[i].equals("."))
                        v.setRsId(data[i]);
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
                                        v.setVariantType("snp");
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
                        return new ArrayList<>(0);
                    }
                    break;
            }
        }
        v.setGenicStatus( isGenic(v) ? "GENIC":"INTERGENIC"  );
        v.setMapKey(mapKey);
        v.setSpeciesTypeKey(3);
        List<VariantMapData> dbVars = dao.getVariant(v);
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
                int varFreq = Integer.parseInt(depths[i+1]);
                if (varFreq==0)
                    continue;
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
                        var = var.replaceFirst(v.getPaddingBase(),"");
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
                        int sampleCnt = dao.getVariantSampleDetailCount((int)dbVar.getId(),s.getId());
                        if (sampleCnt==0){
                            VariantSampleDetail newSample = new VariantSampleDetail();
                            newSample.setId(dbVar.getId());
                            newSample.setSampleId(s.getId());
                            newSample.setDepth(totalDepth);
//                            int varFreq = Integer.parseInt(depths[i+1]);
                            newSample.setVariantFrequency(varFreq);
                            zygosity.computeZygosityStatus(newSample.getVariantFrequency(),newSample.getDepth(),s.getGender(),dbVar, newSample);

                            int zygPercentRead = varFreq / newSample.getDepth();
                            newSample.setZygosityPercentRead(zygPercentRead);
                            samples.add(newSample);

                        }
                        if (dbVar.getEndPos() != copy.getEndPos() && copy.getEndPos()!=0) {
                            dbVar.setEndPos(copy.getEndPos());
                            tobeUpdated.add(dbVar);
                        }
                    }
                }// end check with database vars
                if (!exist) {
                    RgdId r = dao.createRgdId(RgdId.OBJECT_KEY_VARIANTS, "ACTIVE", "created by HRDP Load Pipeline", mapKey);
                    copy.setId(r.getRgdId());
//                    var.setId(r.getRgdId());
                    // create sample
                    VariantSampleDetail newSample = new VariantSampleDetail();
                    newSample.setId(copy.getId());
                    newSample.setSampleId(s.getId());
                    newSample.setDepth(totalDepth);

                    newSample.setVariantFrequency(varFreq);
                    zygosity.computeZygosityStatus(newSample.getVariantFrequency(), newSample.getDepth(), s.getGender(), v, newSample);

                    int zygPercentRead = varFreq / newSample.getDepth();
                    newSample.setZygosityPercentRead(zygPercentRead);
                    samples.add(newSample);
                    newVars.add(copy);

                }

            }

            vars.addAll(newVars);
        }
        else {
            // check if variant exists, if not generate rgdId
            // check if sample exists only if variant exists
            int varFreq = Integer.parseInt(depths[1]);
            if (varFreq!=0) {
                boolean exist = false;
                for (VariantMapData dbVar : dbVars) {
                    if (Utils.stringsAreEqual(v.getReferenceNucleotide(), dbVar.getReferenceNucleotide()) && Utils.stringsAreEqual(v.getVariantNucleotide(), dbVar.getVariantNucleotide())
                            && v.getStartPos() == dbVar.getStartPos()) {
                        exist = true;

                        int sampleCnt = dao.getVariantSampleDetailCount((int) dbVar.getId(), s.getId());
                        if (sampleCnt == 0) {
                            VariantSampleDetail newSample = new VariantSampleDetail();
                            newSample.setId(dbVar.getId());
                            newSample.setSampleId(s.getId());
                            newSample.setDepth(totalDepth);
//                        int varFreq = Integer.parseInt(depths[1]);
                            newSample.setVariantFrequency(varFreq);
                            zygosity.computeZygosityStatus(newSample.getVariantFrequency(), newSample.getDepth(), s.getGender(), dbVar, newSample);

                            int zygPercentRead = varFreq / newSample.getDepth();
                            newSample.setZygosityPercentRead(zygPercentRead);
                            samples.add(newSample);

                        }
                        if (dbVar.getEndPos() != v.getEndPos() && v.getEndPos() != 0) {
                            dbVar.setEndPos(v.getEndPos());
                            tobeUpdated.add(dbVar);
                        }
                    }
                }
                if (!exist) {

                    RgdId r = dao.createRgdId(RgdId.OBJECT_KEY_VARIANTS, "ACTIVE", "created by HRDP Load Pipeline", mapKey);
                    v.setId(r.getRgdId());
                    // create sample
                    VariantSampleDetail newSample = new VariantSampleDetail();
                    newSample.setId(v.getId());
                    newSample.setSampleId(s.getId());
                    newSample.setDepth(totalDepth);
//                    int varFreq = Integer.parseInt(depths[1]);
                    newSample.setVariantFrequency(varFreq);
                    zygosity.computeZygosityStatus(newSample.getVariantFrequency(), newSample.getDepth(), s.getGender(), v, newSample);

                    int zygPercentRead = varFreq / newSample.getDepth();
                    newSample.setZygosityPercentRead(zygPercentRead);
                    samples.add(newSample);
                    vars.add(v);

                }
//            vars.add(v);
            }
        }
        return vars;
    }

    Integer getStrainRgdId(String sampleName) throws Exception {
        int strainStop = sampleName.indexOf('(')-1;
        String strainName = sampleName.substring(0,strainStop);
        return dao.getStrainRgdIdByTaglessStrainSymbol(strainName);
    }

    private BufferedReader openFile(String fileName) throws IOException {

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

    boolean isGenic(VariantMapData v) throws Exception {

        GeneCache geneCache = geneCacheMap.get(v.getChromosome());
        if( geneCache==null ) {
            geneCache = new GeneCache();
            geneCacheMap.put(v.getChromosome(), geneCache);
            geneCache.loadCache(mapKey, v.getChromosome(), DataSourceFactory.getInstance().getDataSource());
        }
        List<Integer> geneRgdIds = geneCache.getGeneRgdIds((int)v.getStartPos(), (int)v.getEndPos());
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
}
