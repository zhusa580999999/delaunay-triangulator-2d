# Microsoft Developer Studio Project File - Name="spdelaunay3d" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=spdelaunay3d - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "spdelaunay3d.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "spdelaunay3d.mak" CFG="spdelaunay3d - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "spdelaunay3d - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "spdelaunay3d - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "spdelaunay3d - Win32 Release"

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
# ADD LINK32 lib/SPlib.lib lib/SVlib.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib /nologo /subsystem:console /machine:I386
# SUBTRACT LINK32 /profile
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Cmds=copy Release\spdelaunay3d.exe spdelaunay3d.exe
# End Special Build Tool

!ELSEIF  "$(CFG)" == "spdelaunay3d - Win32 Debug"

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
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I "inc" /I "stl" /D "_DEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 libD/SPlib.lib libD/SVlib.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# SUBTRACT LINK32 /profile

!ENDIF 

# Begin Target

# Name "spdelaunay3d - Win32 Release"
# Name "spdelaunay3d - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\src\delaunay3.cpp
# End Source File
# Begin Source File

SOURCE=.\src\delaunay3d.cpp
# End Source File
# Begin Source File

SOURCE=.\src\fopengzipped.cpp
# End Source File
# Begin Source File

SOURCE=.\src\predicates.cpp
# End Source File
# Begin Source File

SOURCE=.\src\spdelaunay3d.cpp
# End Source File
# Begin Source File

SOURCE=.\src\spreadscattered.cpp
# End Source File
# Begin Source File

SOURCE=.\src\sscontainer3d.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\inc\delaunay3.h
# End Source File
# Begin Source File

SOURCE=.\inc\spdelaunay3d.h
# End Source File
# Begin Source File

SOURCE=.\src\spdelaunay3d.h
# End Source File
# Begin Source File

SOURCE=.\inc\spreader.h
# End Source File
# Begin Source File

SOURCE=.\inc\spreader_node.h
# End Source File
# Begin Source File

SOURCE=.\inc\spreader_spa.h
# End Source File
# Begin Source File

SOURCE=.\inc\spreader_spb.h
# End Source File
# Begin Source File

SOURCE=.\inc\spreadscattered.h
# End Source File
# Begin Source File

SOURCE=.\inc\spwriter.h
# End Source File
# Begin Source File

SOURCE=.\inc\sscontainer3d.h
# End Source File
# Begin Source File

SOURCE=.\inc\svwriter.h
# End Source File
# Begin Source File

SOURCE=.\inc\svwriter_nil.h
# End Source File
# Begin Source File

SOURCE=.\inc\svwriter_sva.h
# End Source File
# Begin Source File

SOURCE=.\inc\svwriter_svb.h
# End Source File
# Begin Source File

SOURCE=.\inc\vec3fv.h
# End Source File
# Begin Source File

SOURCE=.\inc\vec3iv.h
# End Source File
# Begin Source File

SOURCE=.\inc\vecnfv.h
# End Source File
# Begin Source File

SOURCE=.\inc\vecniv.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
