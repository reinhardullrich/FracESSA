import os

# print(os.getcwd())

print('-----------------------------------------------------------LINUX------------------------------------------------------')
os.system('rm -r /home/reinhard/REF/REF/build_linux64/*')
os.system('cmake -S /home/reinhard/REF/REF/ -B /home/reinhard/REF/REF/build_linux64 -DCMAKE_BUILD_TYPE=Release')
os.chdir(os.path.dirname("/home/reinhard/REF/REF/build_linux64/"))  # change cwd to the directory where the file is!
os.system('make VERBOSE=1')


print('-----------------------------------------------------------WINDOWS------------------------------------------------------')
os.system('rm -r /home/reinhard/REF/REF/build_win64/*')
os.system('cmake --debug -S /home/reinhard/REF/REF/ -B /home/reinhard/REF/REF/build_win64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../win64.cmake')
os.chdir(os.path.dirname("/home/reinhard/REF/REF/build_win64/"))  # change cwd to the directory where the file is!
os.system('make VERBOSE=1')
