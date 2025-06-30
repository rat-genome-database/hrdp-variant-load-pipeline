package edu.mcw.rgd.hrdpVariantLoad;

import edu.mcw.rgd.dao.DataSourceFactory;
import edu.mcw.rgd.datamodel.Gene;
import edu.mcw.rgd.datamodel.MapData;
import edu.mcw.rgd.datamodel.Sample;
import edu.mcw.rgd.datamodel.variants.VariantMapData;
import edu.mcw.rgd.datamodel.variants.VariantSampleDetail;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.*;

public class GenicQc {
    private String version;
    private String inputDir;
    private int mapKey;

    protected Logger logger = LogManager.getLogger("status");
    private DAO dao = new DAO();
    private SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

    void run() throws Exception{
        logger.info(getVersion());
        logger.info("   "+dao.getConnection());

        long pipeStart = System.currentTimeMillis();
        logger.info("Pipeline started at "+sdt.format(new Date(pipeStart))+"\n");

        File folder = new File(inputDir);
        ArrayList<File> files = new ArrayList<>();
        dao.listFilesInFolder(folder, files);
        for (File file : files) {
            parse(file);
        }

        logger.info("Total pipeline runtime -- elapsed time: "+
                Utils.formatElapsedTime(pipeStart,System.currentTimeMillis()));
    }

    void parse(File file) throws Exception{

        // parse file and insert/add samples to variants
        logger.info("\tBegin parsing file... " + file.getName());
        BufferedReader br = dao.openFile(file.getAbsolutePath());
        String lineData;
        List<VariantMapData> tobeUpdated = new ArrayList<>();
        while ((lineData = br.readLine()) != null){
            if (lineData.startsWith("#"))
                continue;


            parseLineData(lineData, tobeUpdated);

        } // end of file read

        if (!tobeUpdated.isEmpty()){
            // update genic status
            logger.info("\t\tUpdating genic status: "+tobeUpdated.size());
            dao.updateVariantMapDataGenicStatus(tobeUpdated);
        }

        logger.info("\tEnd parsing file... " + file.getName());
//        System.out.println(variants.size());
    }

    void parseLineData(String lineData, List<VariantMapData> tobeUpdated) throws Exception {
        VariantMapData v = new VariantMapData();
        List<VariantMapData> vars = new ArrayList<>();
        boolean needCopyVar = false;
        boolean needCopyRef = false;
        Integer totalDepth = null;
        String[] data = lineData.split("\t");
        String[] depths = new String[0];
        // first check if ref or alt has a ','
        // if yes, make a copy of
        for (int i = 0; i < data.length; i++) {
            switch (i) {
                case 0: // chrom
                    // chrM is for MT
                    if (data[i].contains("unplaced") || data[i].contains("unloc") || data[i].contains("contig") || data[i].contains("scaffold")) {
                        return;
                    }
                    v.setChromosome(data[i].replace("chr", ""));
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
                    } else {
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
                    } else {

                        if (data[i].equals("*")) {
                            v.setVariantNucleotide(null);
                            v.setEndPos(v.getStartPos() + v.getReferenceNucleotide().length());
                            v.setVariantType("deletion");
                        } else {
                            String var = data[i];
                            v.setVariantNucleotide(var);
                            if (v.getReferenceNucleotide().length() > var.length() && var.length() == 1) {
                                // deletion
                                v.setPaddingBase(var);
                                v.setStartPos(v.getStartPos() + 1);
                                v.setVariantNucleotide(null);
                                String ref = v.getReferenceNucleotide().substring(1);
                                v.setReferenceNucleotide(ref);
                                v.setEndPos(v.getStartPos() + v.getReferenceNucleotide().length());
                                v.setVariantType("deletion");
                            } else if (v.getReferenceNucleotide().length() > var.length() && v.getReferenceNucleotide().startsWith(var)) {
                                v.setPaddingBase(var);
                                v.setStartPos(v.getStartPos() + var.length());
                                v.setVariantNucleotide(null);
                                String ref = v.getReferenceNucleotide().replaceFirst(var, "");
                                v.setReferenceNucleotide(ref);
                                v.setEndPos(v.getStartPos() + ref.length());
                                v.setVariantType("deletion");
                            } else if (var.length() > v.getReferenceNucleotide().length() && v.getReferenceNucleotide().length() == 1) {
                                // insertion
                                v.setPaddingBase(v.getReferenceNucleotide());
                                v.setStartPos(v.getStartPos() + 1);
                                v.setEndPos(v.getStartPos() + 1);
                                v.setReferenceNucleotide(null);
                                var = var.substring(1);
                                v.setVariantNucleotide(var);
                                v.setVariantType("insertion");
                            } else if (var.length() > v.getReferenceNucleotide().length() && var.startsWith(v.getReferenceNucleotide())) {
                                v.setPaddingBase(v.getReferenceNucleotide());
                                v.setStartPos(v.getStartPos() + v.getReferenceNucleotide().length());
                                v.setEndPos(v.getStartPos() + 1);
                                v.setReferenceNucleotide(null);
                                var = var.replaceFirst(v.getPaddingBase(), "");
                                v.setVariantNucleotide(var);
                                v.setVariantType("insertion");
                            } else {
                                v.setVariantNucleotide(data[i]);
                                if (v.getReferenceNucleotide().length() == v.getVariantNucleotide().length()) {
                                    v.setEndPos(v.getStartPos() + 1);
                                    if (v.getReferenceNucleotide().length() > 1) {
                                        v.setVariantType("mnv");
                                        v.setEndPos(v.getStartPos() + v.getReferenceNucleotide().length());
                                    } else {
                                        v.setVariantType("snp");
                                        v.setEndPos(v.getStartPos() + 1);
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
                    break;
            }
        } // end for
//            if (Utils.stringsAreEqual(v.getRsId(), "rs3319176509"))
//                System.out.println("here");
            List<VariantMapData> variantsInTable = dao.getVariantsWithGeneLocation(mapKey, v.getChromosome(), (int) v.getStartPos(), (int) v.getEndPos());
            for (VariantMapData variant : variantsInTable) {
                String oldGenicStat = variant.getGenicStatus();
                    List<MapData> mapData = dao.getMapDataWithinRange((int)variant.getStartPos(),(int)variant.getEndPos(),variant.getChromosome(),mapKey,1);
                    List<MapData> dataList = new ArrayList<>();
                    if (mapData.size()>0) {
                        dataList.addAll(mapData);
                    }
                    if (!dataList.isEmpty()){
                        variant.setGenicStatus("GENIC");
                    }
                    else {
                        variant.setGenicStatus("INTERGENIC");
                    }
//                if (isGenic(variant)) {
//                    variant.setGenicStatus("GENIC");
//                } else {
//                    variant.setGenicStatus("INTERGENIC");
//                }
                if (!Utils.stringsAreEqualIgnoreCase(variant.getGenicStatus(), oldGenicStat))
                    tobeUpdated.add(variant);
                break;
        } // end var for

    }


    boolean isGenic(VariantMapData v) throws Exception {

        GeneCache geneCache = geneCacheMap.get(v.getChromosome());
        if( geneCache==null ) {
            geneCache = new GeneCache();
            geneCacheMap.put(v.getChromosome(), geneCache);
            geneCache.loadCache(mapKey, v.getChromosome(), DataSourceFactory.getInstance().getDataSource());
        }
        List<Integer> geneRgdIds = geneCache.getGeneRgdIds((int)v.getStartPos(), (int)v.getStartPos());
        if (Utils.stringsAreEqual(v.getRsId(), "rs3319176509")){
            logger.info("rs3319176509 - a problem child that is genic when it is not");
            for (int id : geneRgdIds){
                logger.info(id);
            }
        }
        return !geneRgdIds.isEmpty();
    }
    Map<String, GeneCache> geneCacheMap = new HashMap<>();

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }

    public void setInputDir(String inputDir) {
        this.inputDir = inputDir;
    }

    public String getInputDir() {
        return inputDir;
    }

    public void setMapKey(int mapKey) {
        this.mapKey = mapKey;
    }

    public int getMapKey() {
        return mapKey;
    }
}
