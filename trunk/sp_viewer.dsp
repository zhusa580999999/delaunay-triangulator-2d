# Microsoft Developer Studio Project File - Name="sp_viewer" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=sp_viewer - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "sp_viewer.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "sp_viewer.mak" CFG="sp_viewer - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "sp_viewer - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "sp_viewer - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "sp_viewer - Win32 Release"

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
# ADD CPP /nologo /W3 /GX /O2 /I "inc" /I "stl" /I "src" /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 lib/SPlib.lib lib/SMlib.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Cmds=copy Release\sp_viewer.exe sp_viewer.exe
# End Special Build Tool

!ELSEIF  "$(CFG)" == "sp_viewer - Win32 Debug"

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
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I "inc" /I "stl" /I "src" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D KEY_TYPE=EXTu32 /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 libD/SPlib.lib libD/SMlib.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "sp_viewer - Win32 Release"
# Name "sp_viewer - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\src\delaunay2.cpp
# End Source File
# Begin Source File

SOURCE=.\src\delaunay3.cpp
# End Source File
# Begin Source File

SOURCE=.\src\fopengzipped.cpp
# End Source File
# Begin Source File

SOURCE=.\src\predicates.cpp
# End Source File
# Begin Source File

SOURCE=.\src\sp_viewer.cpp
# End Source File
# Begin Source File

SOURCE=.\src\spdelaunay2d.cpp
# End Source File
# Begin Source File

SOURCE=.\src\spdelaunay3d.cpp
# End Source File
# Begin Source File

SOURCE=.\src\spreadscattered.cpp
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

SOURCE=.\src\Delaunay2.h
# End Source File
# Begin Source File

SOURCE=.\src\delaunay3.h
# End Source File
# Begin Source File

SOURCE=..\src\integercompressor_newer.h
# End Source File
# Begin Source File

SOURCE=..\src\positionquantizer_new.h
# End Source File
# Begin Source File

SOURCE=..\inc\smreader.h
# End Source File
# Begin Source File

SOURCE=..\inc\smreader_ply.h
# End Source File
# Begin Source File

SOURCE=..\inc\smreader_sma.h
# End Source File
# Begin Source File

SOURCE=..\inc\smreader_smb.h
# End Source File
# Begin Source File

SOURCE=..\inc\spconverter.h
# End Source File
# Begin Source File

SOURCE=.\src\spdelaunay2d.h
# End Source File
# Begin Source File

SOURCE=.\src\spdelaunay3d.h
# End Source File
# Begin Source File

SOURCE=..\inc\spreader.h
# End Source File
# Begin Source File

SOURCE=..\inc\spreader_node.h
# End Source File
# Begin Source File

SOURCE=..\inc\spreader_spa.h
# End Source File
# Begin Source File

SOURCE=..\inc\spreader_spb.h
# End Source File
# Begin Source File

SOURCE=..\inc\spreader_spc.h
# End Source File
# Begin Source File

SOURCE=..\inc\spreadscattered.h
# End Source File
# Begin Source File

SOURCE=..\inc\sscontainer2d.h
# End Source File
# Begin Source File

SOURCE=..\inc\sscontainer3d.h
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
