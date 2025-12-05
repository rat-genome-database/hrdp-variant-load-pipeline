package edu.mcw.rgd.hrdpVariantLoad;

import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;


public class Manager {

    public static void main(String[] args) throws Exception {
        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        try {
            for (int i = 0; i < args.length; i++){
                switch (args[i]){
                    case "--runLoad":
                        HrdpVariants hrdpVariants = (HrdpVariants) (bf.getBean("hrdpVariants"));
                        hrdpVariants.main();
                        break;
                    case "--genicQc":
                        GenicQc qc = (GenicQc) (bf.getBean("genicQc"));
                        qc.run();
                        break;
                }
            }

        }
        catch (Exception e) {
            Utils.printStackTrace(e, LogManager.getLogger("status"));
            throw e;
        }
    }

}