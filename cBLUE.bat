ECHO off

SET where=C:\Windows\System32\where.exe
SET anaconda_config=%~dp0%anaconda_dir.txt

:: create anaconda_config file if it doesn't exist
IF NOT EXIST %anaconda_config% (
	ECHO "cBLUE doesn't know where your Anaconda directory is yet."
	ECHO "Locating Anaconda directory...
	ECHO "(result will be stored in anaconda_dir.txt for future use)"
	for /f "delims=_ tokens=1" %%a in ('%where% /r %userprofile% _conda.exe') do ECHO %%a >anaconda_dir.txt
) 

:: read anaconda dir from anaconda_config file
for /f "delims= " %%x in (%anaconda_config%) do set anaconda_dir=%%x
SET conda=%anaconda_dir%condabin\conda.bat
SET env_dir=%anaconda_dir%envs
SET env_name=cblue

:: activate cBLUE conda environment (create conda environment if it doesn't exist)
ECHO "activating cBLUE conda environment..."
CALL %conda% activate %env_name%
IF %ERRORLEVEL% == 0 (
	%env_dir%\%env_name%\python.exe %~dp0%CBlueApp.py
) ELSE (
	ECHO "Oh nooo, a valid cBLUE environment doesn't exist.  Let's create one..."
	%conda% env create --file %~dp0%cBLUE.yml

	ECHO "There, you've now got a cBLUE environment. (You won't have to do this part again...for this particular version)"
	ECHO "activating environment..."
	%conda% activate %env_name%
	%env_dir%\%env_name%\python.exe %~dp0%CBlueApp.py
)

cmd /k