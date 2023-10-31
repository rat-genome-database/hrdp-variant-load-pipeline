package edu.mcw.rgd.hrdpVariantLoad;

import edu.mcw.rgd.dao.impl.SampleDAO;
import edu.mcw.rgd.dao.impl.StrainDAO;
import edu.mcw.rgd.dao.impl.variants.VariantDAO;
import edu.mcw.rgd.datamodel.Sample;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

public class HrdpVariants {

    private String version;
    int sampleIdStart;
    private String vcfPath = null;
    int mapKey = 0;

    public void main(String[] args) throws Exception {

        for (int i = 0; i < args.length; i++){
            switch (args[i]){
                case "--inputDir":
                    vcfPath = args[++i];
                    break;
                case "--mapKey":
                    mapKey = Integer.parseInt(args[++i]);
                    break;

            }
        }

        // loops through files
        File folder = new File(vcfPath);
        ArrayList<File> files = new ArrayList<>();
        listFilesInFolder(folder, files);

        for (File file : files) {
            parse(file);
        }
    }

    void parse(File file) throws Exception{
        // create samples and insert them
        // parse through files and insert variants
        String strain = getStrainName(file.getName());
        System.out.println(strain);
        Integer strainRgdId = getStrainRgdId(strain);

        SampleDAO sampleDAO = new SampleDAO();
        Sample sample = sampleDAO.getSampleByStrainRgdIdNMapKey(strainRgdId,mapKey);

        if (sample == null) {
            sample = new Sample();
            sample.setId(sampleIdStart);
            sampleIdStart++;
            sample.setAnalysisName(strain);

            if (strainRgdId != 0) {
                sample.setStrainRgdId(strainRgdId);
            }
            sample.setGender("U");
            sample.setDescription("Dr. Mindy Dwinell - Hybrid rat diversity program");
            sample.setPatientId(372); // default is rat 7.2
            sample.setMapKey(mapKey);
            sample.setGrantNumber("R24OD022617");
            VariantDAO vdao = new VariantDAO();
            vdao.insertSample(sample);
        }

        // parse file and insert/add samples to variants

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

    Integer getStrainRgdId(String sampleName) throws Exception {
        StrainDAO sdao = new StrainDAO();
        int strainStop = sampleName.indexOf('(')-1;
        String strainName = sampleName.substring(0,strainStop);
        return sdao.getStrainRgdIdByTaglessStrainSymbol(strainName);
    }
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
}
