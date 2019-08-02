public class TestProgram{
	 // Test Driver
   public static void main(String[] args) {
   		F5C f5c = new F5C();

		int result = f5c.init("f5c index -d ../test/ecoli_2kb_region/fast5_files/ ../test/ecoli_2kb_region/reads.fasta");  // Create an instance and invoke the native method
		System.out.println("ran test program. result: "+result);  
	}
}