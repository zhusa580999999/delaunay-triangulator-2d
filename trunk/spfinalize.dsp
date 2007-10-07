# Microsoft Developer Studio Project File - Name="spfinalizer" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=spfinalizer - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "spfinalize.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "spfinalize.mak" CFG="spfinalizer - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "spfinalizer - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "spfinalizer - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "spfinalizer - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "inc" /I "stl" /D "NDEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 lib/SPlib.lib lib/SMlib.lib lib/SVlib.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Cmds=copy Release\spfinalize.exe spfinalize.exe
# End Special Build Tool

!ELSEIF  "$(CFG)" == "spfinalizer - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I "inc" /I "stl" /D "_DEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /FR /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 libD/SPlib.lib libD/SMlib.lib libD/SVlib.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "spfinalizer - Win32 Release"
# Name "spfinalizer - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\src\fopengzipped.cpp
# End Source File
# Begin Source File

SOURCE=.\src\spfinalizer.cpp
# End Source File
# Begin Source File

SOURCE=.\src\sscontainer2d.cpp
# End Source File
# Begin Source File

SOURCE=.\src\sscontainer3d.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\src\positionquantizer_new.h
# End Source File
# Begin Source File

SOURCE=..\inc\smreader.h
# End Source File
# Begin Source File

SOURCE=..\inc\smreader_sma.h
# End Source File
# Begin Source File

SOURCE=..\inc\smreader_smb.h
# End Source File
# Begin Source File

SOURCE=..\inc\smreader_smc.h
# End Source File
# Begin Source File

SOURCE=..\inc\spconverter.h
# End Source File
# Begin Source File

SOURCE=..\inc\spreader.h
# End Source File
# Begin Source File

SOURCE=..\inc\spreader_node.h
# End Source File
# Begin Source File

SOURCE=..\inc\spreader_ply.h
# End Source File
# Begin Source File

SOURCE=..\inc\spreader_raw.h
# End Source File
# Begin Source File

SOURCE=..\inc\spreader_raw_d.h
# End Source File
# Begin Source File

SOURCE=..\inc\spreader_spa.h
# End Source File
# Begin Source File

SOURCE=..\inc\spreader_spb.h
# End Source File
# Begin Source File

SOURCE=..\inc\spreader_tiles.h
# End Source File
# Begin Source File

SOURCE=..\inc\spwriter.h
# End Source File
# Begin Source File

SOURCE=..\inc\spwriter_nil.h
# End Source File
# Begin Source File

SOURCE=..\inc\spwriter_spa.h
# End Source File
# Begin Source File

SOURCE=..\inc\spwriter_spb.h
# End Source File
# Begin Source File

SOURCE=..\inc\sscontainer2d.h
# End Source File
# Begin Source File

SOURCE=..\inc\sscontainer3d.h
# End Source File
# Begin Source File

SOURCE=..\inc\svreader.h
# End Source File
# Begin Source File

SOURCE=..\inc\svreader_sva.h
# End Source File
# Begin Source File

SOURCE=..\inc\svreader_svb.h
# End Source File
# Begin Source File

SOURCE=..\inc\svreader_svc.h
# End Source File
# Begin Source File

SOURCE=..\inc\vec3dv.h
# End Source File
# Begin Source File

SOURCE=..\inc\vec3fv.h
# End Source File
# Begin Source File

SOURCE=..\inc\vec3iv.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
