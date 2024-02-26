import platform
import os
from pathlib import Path

os.chdir(os.path.dirname(os.path.realpath(__file__)))

# if platform.system() == "Windows":

#     # test if local windows ins installed

#     if not os.system('"jre\\javawindows\\bin\\java -version"'):
#         JAVA_AVAILABLE = True
#         WINDOWS_OS = True
#         JAVA_COMMAND = '"jre\\javawindows\\bin\\java -classpath predictorc.jar;cdk-2.7.1.jar;. NewTest mol.mol > mol.csv"'
#         print("\nWINDOWS Local JAVA is available\n")
#     else:
#         JAVA_AVAILABLE = False
#         JAVA_COMMAND = ""
#         print("JAVA is not available")
# elif platform.system() == "Linux":
#     LINUX_OS = True
#     if not os.system('"jre\\amazon-corretto-17.0.4.9.1-linux-x64\\bin\\java -version"'):
#         JAVA_AVAILABLE = True
#         JAVA_COMMAND = '"amazon-corretto-17.0.4.9.1-linux-x64\\bin\\java -classpath predictorc.jar:cdk-2.7.1.jar:. NewTest mol.mol > mol.csv"'
#         print("Linux Local JAVA is available")
#     else:
#         JAVA_AVAILABLE = False
#         JAVA_COMMAND = ""
#         print("JAVA is not available")
# elif platform.system() == "Darwin":
#     MAC_OS = True
#     if not os.system("jre/amazon-corretto-17.jdk/Contents/Home/bin/java --version"):
#         JAVA_AVAILABLE = True
#         JAVA_COMMAND = "jre/amazon-corretto-17.jdk/Contents/Home/bin/java -classpath predictorc.jar:cdk-2.7.1.jar:. NewTest mol.mol > mol.csv"
#         print("MAC Local JAVA is available")
#     else:
#         JAVA_AVAILABLE = False
#         JAVA_COMMAND = ""
#         print("JAVA is not available")

def find_java_commands(wherewhich_cmd):
    java_exes = os.popen(wherewhich_cmd).read().split("\n")
    possible_java_cmds = []
    for cmd in java_exes:
        if cmd == "":
            continue
        pth = Path(cmd)
        qpth = f'"{pth}"'
        print(pth.is_file())
        # print(cmd, os.system(f"{qpth} --version"))
        txt = os.popen(f"{qpth} --version").read()
        # if txt:
            # print(cmd, "txt:", len(txt))
        possible_java_cmds.append(f"{qpth} -classpath predictorc.jar;cdk-2.7.1.jar;. NewTest mol.mol > mol.csv")

        return possible_java_cmds


if platform.system() == "Windows":

    # test if local windows ins installed

    if not os.system('"jre\\javawindows\\bin\\java -version"'):
        JAVA_AVAILABLE = True
        WINDOWS_OS = True
        JAVA_COMMAND = '"jre\\javawindows\\bin\\java -classpath predictorc.jar;cdk-2.7.1.jar;. NewTest mol.mol > mol.csv"'
        print("\nWINDOWS Local JAVA is available\n")
    else:
        possible_java_cmds = find_java_commands("where java")
        if len(possible_java_cmds) == 0:
            JAVA_AVAILABLE = False
            JAVA_COMMAND = ""
            print("JAVA is not available")
        else:
            JAVA_AVAILABLE = True
            WINDOWS_OS = True
            JAVA_COMMAND = possible_java_cmds[-1]
            print(JAVA_COMMAND)
            print("WINDOWS System JAVA is available")

elif platform.system() == "Linux":
    LINUX_OS = True
    if not os.system('"jre\\amazon-corretto-17.0.4.9.1-linux-x64\\bin\\java -version"'):
        JAVA_AVAILABLE = True
        JAVA_COMMAND = '"amazon-corretto-17.0.4.9.1-linux-x64\\bin\\java -classpath predictorc.jar:cdk-2.7.1.jar:. NewTest mol.mol > mol.csv"'
        print("Linux Local JAVA is available")
    else:
        possible_java_cmds = find_java_commands("which java")
        if len(possible_java_cmds) == 0:
            JAVA_AVAILABLE = False
            JAVA_COMMAND = ""
            print("JAVA is not available")
        else:
            JAVA_AVAILABLE = True
            WINDOWS_OS = True
            JAVA_COMMAND = possible_java_cmds[-1]
            print("Linux System JAVA is available")
elif platform.system() == "Darwin":
    MAC_OS = True
    if not os.system("jre/amazon-corretto-17.jdk/Contents/Home/bin/java --version"):
        JAVA_AVAILABLE = True
        JAVA_COMMAND = "jre/amazon-corretto-17.jdk/Contents/Home/bin/java -classpath predictorc.jar:cdk-2.7.1.jar:. NewTest mol.mol > mol.csv"
        print("MAC Local JAVA is available")
    else:
        possible_java_cmds = find_java_commands("which java")
        if len(possible_java_cmds) == 0:
            JAVA_AVAILABLE = False
            JAVA_COMMAND = ""
            print("JAVA is not available")
        else:
            JAVA_AVAILABLE = True
            WINDOWS_OS = True
            JAVA_COMMAND = possible_java_cmds[-1]
            print("MacOs System JAVA is available")
