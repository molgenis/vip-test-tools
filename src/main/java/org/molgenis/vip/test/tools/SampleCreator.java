package org.molgenis.vip.test.tools;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class SampleCreator {

  private static final String PREFIX = "SAMPLE";
  private Random random = new Random();

  public void run(String[] args) {
    CommandLineParser commandLineParser = new DefaultParser();

    Options options = new Options();
    options.addOption(
        Option.builder("i").longOpt("input").hasArg(true).desc("").required().build());
    options.addOption(Option.builder("o").longOpt("output").hasArg(true).required().desc("").build());
    options.addOption(Option.builder("f").longOpt("force").desc("").build());
    options.addOption(Option.builder("sa").longOpt("samples").desc("nr of samples").hasArg(true).required().build());
    options.addOption(Option.builder("ph").longOpt("phases").desc("yes, no, mixed").hasArg(true).required().build());
    options.addOption(Option.builder("fr").longOpt("frequency").desc("Alternative allele frequency").hasArg(true).required().build());

    CommandLine commandLine;
    try {
      commandLine = commandLineParser.parse(options, args);
    } catch (ParseException exp) {
      System.err.println(exp.getMessage());
      return;
    }
    String inputPath = commandLine.getOptionValue("i");
    File inputFile = new File(inputPath);
    if (!inputFile.exists()) {
      System.err.println(inputFile + " does not exist");
    }
    File outputFile;
    if (commandLine.hasOption("o")) {
      outputFile = new File(commandLine.getOptionValue("o"));
      if (outputFile.exists()) {
        if (!commandLine.hasOption("f")) {
          System.out.println("Output file already exists, use '-f' to overwrite.");
        }
      }
    }else{
      outputFile = new File(inputPath.replace(".vcf",".out.vcf"));
    }
    int nrOfSamples = Integer.parseInt(commandLine.getOptionValue("sa"));
    PhasedModeEnum phasedMode = toPhasedMode(commandLine.getOptionValue("ph"));
    int altFreqency = Integer.parseInt(commandLine.getOptionValue("fr"));

    VCFFileReader reader = new VCFFileReader(inputFile, false);
    List<String> sampleNames = createNames(nrOfSamples);
    VCFHeader header = new VCFHeader(reader.getFileHeader().getMetaDataInInputOrder(), sampleNames);

    VariantContextWriterBuilder vcWriterBuilder =
        new VariantContextWriterBuilder().clearOptions().setOutputFile(outputFile);
    VariantContextWriter writer = vcWriterBuilder.build();
    writer.writeHeader(header);
    for (CloseableIterator<VariantContext> it = reader.iterator(); it.hasNext(); ) {
      VariantContext context = it.next();
      writer.add(addSamples(context, sampleNames, phasedMode, altFreqency));
    }
    writer.close();
  }

  private List<String> createNames(int nrOfSamples) {
    List<String> names = new ArrayList<>();
    for (int i = 0; i < nrOfSamples; i++){
      names.add(PREFIX +i);
    }
    return names;
  }

  private VariantContext addSamples(VariantContext context, List<String> sampleNames,
      PhasedModeEnum phasedMode, int altFreqency) {
    VariantContextBuilder builder = new VariantContextBuilder(context);
    List<Genotype> genotypes = new ArrayList<>();
    int nrOfAlts = context.getAlternateAlleles().size();
    for (String sampleName : sampleNames){
      Allele allele1 = getAllele(context, altFreqency, nrOfAlts);
      Allele allele2 = getAllele(context, altFreqency, nrOfAlts);
      GenotypeBuilder genotypeBuilder = new GenotypeBuilder();
      genotypeBuilder.name(sampleName);
      genotypeBuilder.alleles(Arrays.asList(allele1, allele2));
      if(phasedMode == PhasedModeEnum.PHASED){
        genotypeBuilder.phased(true);
      }else if(phasedMode == PhasedModeEnum.NOT_PHASED){
        genotypeBuilder.phased(false);
      }else{
        Random random = new Random();
        if(random.nextBoolean()){
          genotypeBuilder.phased(true);
        }else{
          genotypeBuilder.phased(false);
        }
      }
      genotypes.add(genotypeBuilder.make());
    }
    return builder.genotypes(genotypes).make();
  }

  private Allele getAllele(VariantContext context, int altFreqency, int nrOfAlts) {
    int randomNumber = random.nextInt(100)+1;
    Allele allele;
    if(randomNumber >= altFreqency){
      allele = context.getReference();
    }else{
      if (nrOfAlts > 1) {
        int randomAlleleNumber = random.nextInt(nrOfAlts - 1);
        allele = context.getAlternateAllele(randomAlleleNumber);
      }else{
        allele = context.getAlternateAllele(0);
      }
    }
    return allele;
  }

  private PhasedModeEnum toPhasedMode(String p) {
    switch (p){
      case "yes":
        return PhasedModeEnum.PHASED;
      case "mixed":
        return PhasedModeEnum.MIXED;
      default:
        return PhasedModeEnum.NOT_PHASED;
    }
  }
}
