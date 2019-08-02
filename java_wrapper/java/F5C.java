public class F5C{
	// Used to load the 'native-lib' library on application startup.
    static {
        System.loadLibrary("F5CJava");
    }

    public static native int init(String command);
}