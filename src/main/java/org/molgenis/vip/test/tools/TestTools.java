package org.molgenis.vip.test.tools;

import org.apache.commons.cli.ParseException;

public class TestTools {
  public static void main(String[] args) throws ParseException {
    if(args.length >= 1){
      String application = args[0];
      switch (application){
        case "SampleCreator":
          new SampleCreator().run(args);
          break;
        case "SampleStripper":
          new SampleStripper().run(args);
          break;
        case "BamCreator":
          new BamCreator().run(args);
          break;
        default:
          throw new UnsupportedOperationException(String.format("No such tool: %s", application));
      }
    }else{
      System.out.println("Usage: java TestTools toolname [tool options]");
    }
  }
}
