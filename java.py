
import platform
import os

print(" in java.py")
print("Platform: ", platform.system())
os.chdir(os.path.dirname(os.path.realpath(__file__)))

if platform.system() == "Windows":
    JAVA_AVAILABLE = True
    # test if local windows ins installed
    print('"jre\\javawindows\\bin\\java"', os.path.exists("jre\\javawindows\\bin"))
    print('"jre\\javawindows\\bin\\java"', os.path.exists(r"jre\\javawindows\\bin"))
    print('"jre\\javawindows\\bin\\java"', os.path.exists(r"jre\javawindows\bin"))
    print('"jre\\javawindows\\bin\\java"', os.path.exists("jre\javawindows\bin"))
    print('"jre\\javawindows\\bin\\java"', os.path.exists(r"jre/javawindows/bin"))
    print('"jre\\javawindows\\bin\\java"', os.path.exists("jre/javawindows/bin"))
    if not os.system('"jre\\javawindows\\bin\\java -version"'):

        WINDOWS_OS = True
        JAVA_COMMAND = '"jre\\javawindows\\bin\\java -classpath predictorc.jar;cdk-2.7.1.jar;. NewTest mol.mol > mol.csv"'
        print("\nWINDOWS Local JAVA is available\n")
    else:
        JAVA_AVAILABLE = False
        JAVA_COMMAND = ""
        print("JAVA is not available")
elif platform.system() == "Linux":
    LINUX_OS = True
    if not os.system('"jre\\javalinux\\bin\\java -version"'):
        JAVA_AVAILABLE = True
        JAVA_COMMAND = '"javalinux\\bin\\java -classpath predictorc.jar:cdk-2.7.1.jar:. NewTest mol.mol > mol.csv"'
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