package org.molgenis.vip.test.tools;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class SampleStripper {
  public void run(String[] args) {
    CommandLineParser commandLineParser = new DefaultParser();

    Options options = new Options();
    options.addOption(
        Option.builder("i").longOpt("input").hasArg(true).desc("").required().build());
    options.addOption(Option.builder("o").longOpt("output").hasArg(true).desc("").required().build());
    options.addOption(Option.builder("f").longOpt("force").desc("").build());

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

    VCFFileReader reader = new VCFFileReader(inputFile, false);
    VCFHeader header = new VCFHeader(reader.getFileHeader().getMetaDataInInputOrder(), Collections.emptySet());
    VariantContextWriterBuilder vcWriterBuilder =
        new VariantContextWriterBuilder().clearOptions().setOutputFile(outputFile);
    VariantContextWriter writer = vcWriterBuilder.build();

    writer.writeHeader(header);
    for (CloseableIterator<VariantContext> it = reader.iterator(); it.hasNext(); ) {
      VariantContext context = it.next();
      List<String> attrs = new ArrayList<>();
      attrs.addAll(context.getAttributes().keySet());
      writer.add(new VariantContextBuilder(context).genotypes(Collections.emptyList()).make());
    }
    writer.close();
  }
}
