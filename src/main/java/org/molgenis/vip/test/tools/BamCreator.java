package org.molgenis.vip.test.tools;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CigarUtil;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;
import java.util.Random;
import java.util.UUID;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class BamCreator {

  private final Random random = new Random();
  private static final int DEFAULT_MAX_READS = 100;
  private static final int DEFAULT_FRAGMENT_SIZE = 100;

  public void run(String args[]) {
    CommandLineParser commandLineParser = new DefaultParser();

    Options options = new Options();
    options.addOption(
        Option.builder("o").longOpt("output").hasArg(true).desc("").required().build());
    options.addOption(
        Option.builder("i").longOpt("input").hasArg(true).desc("").required().build());
    options.addOption(
        Option.builder("fa").longOpt("fasta").hasArg(true).desc("").required().build());
    options.addOption(
        Option.builder("f").longOpt("force").desc("").build());
    options.addOption(
        Option.builder("s").longOpt("sample").hasArg(true).desc("").required().build());

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

    File fastaFile = new File(commandLine.getOptionValue("fasta"));
    String sample = commandLine.getOptionValue("sample");
    File outputFile = new File(commandLine.getOptionValue("output"));

    createBam(inputFile, fastaFile, sample, outputFile);
  }

  private void createBam(File inputFile, File fastaFile,
      String sample, File outputFile) {

    IndexedFastaSequenceFile fasta = getIndexedFastaSequenceFile(fastaFile);

    ReferenceSequence nextsequence = fasta.nextSequence();

    SAMFileHeader header = createSamFileHeader(fasta, nextsequence);
    SAMFileWriter bamFileWriter = createSamFileWriter(fastaFile, outputFile, header);
    VCFFileReader reader = new VCFFileReader(inputFile, false);

    for (CloseableIterator<VariantContext> it = reader.iterator(); it.hasNext(); ) {
      VariantContext context = it.next();
      Genotype genotype = context.getGenotype(sample);
      int depth = genotype.getDP();
      if(depth == -1){
        depth = random.nextInt(DEFAULT_MAX_READS);
      }
      int nrOfReads = depth;
      for (int i = 0; i < nrOfReads; i++) {
        createSingleRead(fasta, header, bamFileWriter, context, genotype);
      }
    }
    bamFileWriter.close();
  }

  private void createSingleRead(IndexedFastaSequenceFile fasta, SAMFileHeader header,
      SAMFileWriter bamFileWriter, VariantContext context, Genotype genotype) {
    int size = DEFAULT_FRAGMENT_SIZE;
    String alt = context.getAlternateAllele(0).getBaseString();

    if (alt.length() <= size) {

      int variantPosition = getRandomInt(1, size);
      if (((size - variantPosition) + alt.length()) > size) {
        variantPosition = variantPosition - (alt.length() - 1);
      }

      int startPos = context.getStart() - variantPosition;
      boolean isVariant = false; // HomRef
      if (genotype.isHet()) {
        isVariant = getRandomBoolean(50);
      } else if (genotype.isHomVar()) {
        isVariant = getRandomBoolean(90);
      }

      byte[] bases =
          getBases(context.getContig(), startPos, alt, variantPosition, size, fasta, isVariant);

      String cigar = isVariant ? getCigar(size, variantPosition) : getCigar(size, -1);

      bamFileWriter.addAlignment(
          generateRecord(
              UUID.randomUUID().toString(),
              context.getContig(),
              startPos,
              bases,
              getStrand(),
              getQualities(size),
              header,
              cigar));
    }
  }

  private SAMFileWriter createSamFileWriter(File fastaFile, File outputFile,
      SAMFileHeader header) {
    SAMFileWriterFactory samFileWriterFactory = new SAMFileWriterFactory();
    samFileWriterFactory.setCreateIndex(true);
    SAMFileWriter bamFileWriter = samFileWriterFactory.makeWriter(header, false, outputFile,
        fastaFile);
    return bamFileWriter;
  }

  private SAMFileHeader createSamFileHeader(IndexedFastaSequenceFile fasta,
      ReferenceSequence nextsequence) {
    SAMSequenceDictionary dict = new SAMSequenceDictionary();
    while (nextsequence != null) {
      SAMSequenceRecord sequenceRecord =
          new SAMSequenceRecord(nextsequence.getName(), nextsequence.getBaseString().length());
      dict.addSequence(sequenceRecord);
      nextsequence = fasta.nextSequence();
    }
    SAMFileHeader header = new SAMFileHeader(dict);
    header.setSortOrder(SortOrder.coordinate);
    return header;
  }

  private IndexedFastaSequenceFile getIndexedFastaSequenceFile(File fastaFile) {
    IndexedFastaSequenceFile fasta = null;
    try {
      fasta = new IndexedFastaSequenceFile(fastaFile);
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    }
    return fasta;
  }

  private SAMRecord generateRecord(
      String name,
      String contig,
      int start,
      byte[] bases,
      boolean isNegativeStrand,
      byte[] qualities,
      SAMFileHeader header,
      String cigar) {
    SAMRecord record = new SAMRecord(header);
    record.setReadName(name);
    record.setReferenceName(contig);
    record.setAlignmentStart(start);
    record.setReadBases(bases);
    record.setCigarString(cigar);
    record.setReadNegativeStrandFlag(isNegativeStrand);
    record.setBaseQualities(qualities);
    return record;
  }

  private String getCigar(int size, int variantOffset) {
    //FIXME: insertions and deletions
    char[] value = new char[size];
    for (int i = 0; i < size; i++) {
      if (i == variantOffset) {
        value[i] = 'M';
      } else {
        value[i] = '=';
      }
    }
    return CigarUtil.cigarStringFromArray(value);
  }

  private byte[] getQualities(int size) {
    byte[] value = new byte[size];
    for (int i = 0; i < size; i++) {
      value[i] = (byte) random.nextInt(50);
    }
    return value;
  }

  private int getRandomInt(int from, int to) {
    return from + new Random().nextInt(Math.abs(to - from));
  }

  private byte[] getBases(
      String contig,
      int start,
      String alt,
      int variantPosition,
      int size,
      ReferenceSequenceFile fasta,
      boolean isVariant) {
    byte[] bases = fasta.getSubsequenceAt(contig, start, start + size - 1).getBases();

    if (isVariant) {
      for (int i = 0; i < alt.length(); i++) {
        if (List.of('A', 'T', 'G', 'C').contains(alt.charAt(i))) {
          bases[(variantPosition) + i] = (byte) alt.charAt(i);
        } else {
          bases[(variantPosition) + i] = (byte) 'N';
        }
      }
    }
    return bases;
  }

  private boolean getStrand() {
    return random.nextBoolean();
  }

  public boolean getRandomBoolean(int chance){
    return random.nextInt(100) < chance;
  }
}
