README.txt

Get_seqs_from_fastq

The program get_seqs_from_fastq is compiled with the following five *.cpp and *.h files:
get_seqs_from_fastq.cpp
get_seqs_from_fastq.h
hp_dynamicintegerarra.cppy
hp_dynamicintegerarray.h
hp_dynamicstringarray.cpp
hp_dynamicstringarray.h
hp_quicksort.cpp
hp_quicksort.h
hp_readtextfile.cpp
hp_readtextfile.h

We used versio 4.7.3 of the MingW toolchain (GCC Assembler, GCC C++ compiler, MingW linker)
The files were compiled and linked with the following commands:
g++ -O0 -g3 -Wall -c -fmessage-length=0 -o "src\\hp_dynamicintegerarray.o" "..\\src\\hp_dynamicintegerarray.cpp" 
g++ -O0 -g3 -Wall -c -fmessage-length=0 -o "src\\get_seqs_from_fastq.o" "..\\src\\get_seqs_from_fastq.cpp" 
g++ -O0 -g3 -Wall -c -fmessage-length=0 -o "src\\hp_dynamicstringarray.o" "..\\src\\hp_dynamicstringarray.cpp" 
g++ -O0 -g3 -Wall -c -fmessage-length=0 -o "src\\hp_quicksort.o" "..\\src\\hp_quicksort.cpp" 
g++ -O0 -g3 -Wall -c -fmessage-length=0 -o "src\\hp_readtextfile.o" "..\\src\\hp_readtextfile.cpp" 
g++ -o get_seqs_from_fastq.exe "src\\hp_readtextfile.o" "src\\hp_quicksort.o" "src\\hp_dynamicstringarray.o" "src\\hp_dynamicintegerarray.o" "src\\get_seqs_from_fastq.o"

You must change the string values of "file_in" and "outfilename" in "get_seqs_from_fastq.cpp" to point to your FASTQ file and to a filename and directory where you want to store the output


