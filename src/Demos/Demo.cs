using System;
using ConsensusCore;

public class Demo
{
   public static void Main()
   {
       var vs = ConsensusCore.Version.VersionString();
       System.Console.WriteLine(vs);
   }
}
