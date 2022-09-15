@ECHO OFF
REM .-- Prepare the Command Processor
SETLOCAL ENABLEEXTENSIONS
SETLOCAL ENABLEDELAYEDEXPANSION
 
ECHO.   
ECHO.  Clean ================================================= 
ECHO.   
:menuLOOP
for /f "tokens=1,2,* delims=_ " %%A in ('"findstr /b /c:":menu_" "%~f0""') do echo.  %%B  %%C
set choice=
echo.&set /p choice=Make a choice or hit ENTER to quit: ||GOTO:EOF
echo.&call:menu_%choice%
GOTO:menuLOOP

::-----------------------------------------------------------
:: menu functions follow below here
::-----------------------------------------------------------

:menu_1   delete *.o files
del *.o
GOTO:EOF
 
:menu_2   delete *.mod files
del *.mod
GOTO:EOF

:menu_3 delete *.dat files (except input.dat,license.dat,group.dat,d_electron.dat,f_electron.dat,CCF.dat)
md save
xcopy input.dat save >NUL
xcopy license.dat save >NUL
xcopy group.dat save >NUL
xcopy d_electron.dat save >NUL
xcopy f_electron.dat save >NUL
xcopy CCF.dat save >NUL
del *.dat
move save\*.* .  >NUL
rd save
GOTO:EOF

:menu_4 dir
dir
GOTO:EOF

:menu_  
 
:menu_C   Clear Screen
cls
GOTO:EOF
