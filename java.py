
import platform
import os

os.chdir(os.path.dirname(os.path.realpath(__file__)))

if platform.system() == "Windows":
    
    # test if local windows ins installed

    if not os.system('"jre\\javawindows\\bin\\java -version"'):
        JAVA_AVAILABLE = True
        WINDOWS_OS = True
        JAVA_COMMAND = '"jre\\javawindows\\bin\\java -classpath predictorc.jar;cdk-2.7.1.jar;. NewTest mol.mol > mol.csv"'
        print("\nWINDOWS Local JAVA is available\n")
    else:
        JAVA_AVAILABLE = False
        JAVA_COMMAND = ""
        print("JAVA is not available")
elif platform.system() == "Linux":
    LINUX_OS = True
    if not os.system('"jre\\amazon-corretto-17.0.4.9.1-linux-x64\\bin\\java -version"'):
        JAVA_AVAILABLE = True
        JAVA_COMMAND = '"amazon-corretto-17.0.4.9.1-linux-x64\\bin\\java -classpath predictorc.jar:cdk-2.7.1.jar:. NewTest mol.mol > mol.csv"'
        print("Linux Local JAVA is available")
    else:
        JAVA_AVAILABLE = False
        JAVA_COMMAND = ""
        print("JAVA is not available")
elif platform.system() == "Darwin":
    MAC_OS = True
    if not os.system("jre/amazon-corretto-17.jdk/Contents/Home/bin/java --version"):
        JAVA_AVAILABLE = True
        JAVA_COMMAND = "jre/amazon-corretto-17.jdk/Contents/Home/bin/java -classpath predictorc.jar:cdk-2.7.1.jar:. NewTest mol.mol > mol.csv"
        print("MAC Local JAVA is available")
    else:
        JAVA_AVAILABLE = False
        JAVA_COMMAND = ""
        print("JAVA is not available")