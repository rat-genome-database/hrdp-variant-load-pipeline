package edu.mcw.rgd.hrdpVariantLoad;

import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;

import java.util.List;

public class Manager {

    public static void main(String[] args) throws Exception {
        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        try {
            HrdpVariants hrdpVariants = (HrdpVariants) (bf.getBean("hrdpVariants"));
            hrdpVariants.main(args);
        }
        catch (Exception e) {
            Utils.printStackTrace(e, LogManager.getLogger("status"));
            throw e;
        }
    }

}