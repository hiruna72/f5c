## F5C Java Wrapper 


See [this example](https://github.com/mkowsiak/jnicookbook/blob/master/recipeNo023/Makefile) to understand how to wrap a native library using [JNI](https://www3.ntu.edu.sg/home/ehchua/programming/java/JavaNativeInterface.html)

### Directory structure
```
- java
  |- c
  |- java
  |- lib
  |- target
```        

* __c__ contains the JNI implementation files 
* __java__ contains the java interface for F5C and a TestDriver
* __lib__ contains libf5cshared.so and the [adapter shared lib](http://jnicookbook.owsiak.org/recipe-No-023/)
* __target__ contains executables

### Steps

* Build the adapter shared lib (`libF5CJava.so`)
	* cd to wrappers/java
	*  `g++ -shared -fPIC -I${JAVA_HOME}/include -I${JAVA_HOME}/include/linux -I$PWD/c/include/ c/F5C.cpp lib/libf5cshared.so -o lib/libF5CJava.so`
* Create java Bytecodes:
	* `javac -d target java/F5C.java`
	* `javac -d target java/TestProgram.java -cp "target"`
* Run the test java program
	* `java -d64 -Djava.library.path="$PWD/lib" -cp "target" TestProgram`
* 	If an error is raised and unable to link the library libtwentyface.so export `LD_LIBRARY_PATH` to the appropriate directory. In this case it is `export LD_LIBRARY_PATH=./lib` 

