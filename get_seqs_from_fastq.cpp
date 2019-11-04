#include "get_seqs_from_fastq.h"
using namespace std;

int main()
{

// Open fastq file

char* infile_filepath = new char[1024];
memset(infile_filepath,0,1024);

// Set the following string to the full path to the FASTQ format file from which you want to extract valid reads


strcpy(infile_filepath,"C:\\Users\\hpatterton\\Documents\\Students\\Mzwanele Ngubo\\PhD\\Sequence data\\Sample_H1\\H1_ATCACG_L001_R1_001.fastq");
char* filetype_point = strrchr(infile_filepath,'.');
int number_of_characters = filetype_point-infile_filepath;
char* outfile_filepath = new char[1024];
memset(outfile_filepath,0,1024);
strncpy(outfile_filepath,infile_filepath,number_of_characters);
strcat(outfile_filepath,"_barcodes.txt");

cout << "Fastq file: " <<  infile_filepath << endl;
cout << "Retrieving sequences..." << endl;
// Read all the sequences to an array
HP_DynamicStringArray* string_array = new HP_DynamicStringArray();
HP_ReadTextFile* readfile = new HP_ReadTextFile;

int entries = readfile->ReadFastqFile(infile_filepath, string_array);
cout << entries << " sequences retrieved from fastq file" << endl;
char** primers = new char*[8];
for(int x = 0; x < 8; x++)
{
	primers[x] = new char[18];
	memset(primers[x], 0, 18);
}

/*
 * The primer sequences were copied from Dai et al. (2008) Probing nucleosome function: A highly versatile library of synthetic histone H3 and H4 mutants
 * and correspond to the following:
 *
 * primers[0]: U1h
 * primers[4]: U2h complement
 * primers[1]: U2h
 * primers[5]: U1h complement
 * primers[2]: D1h
 * primers[6]: D2h complement
 * primers[3]: D2h
 * primers[7]: D1h complement
 *
 */

strcpy(primers[0],"ATGTCCACGAGGTCTCT");
strcpy(primers[4], "TACGCTGCAGGTCGAGG");

strcpy(primers[1], "CCTCGACCTGCAGCGTA");
strcpy(primers[5], "AGAGACCTCGTGGACAT");

strcpy(primers[2], "CGGTGTCGGTCTCGTAG");
strcpy(primers[6], "GATGAATTCGAGCTGGG");

strcpy(primers[3], "CCCAGCTCGAATTCATC");
strcpy(primers[7], "CTACGAGACCGACACCG");


int two_primers = 0;
int correct_length = 0;
int one_primer = 0;
int incorrect_length = 0;
int no_primers = 0;

// See if we can match a primer

int primer_length = 17;
int insert_length = 20;
char* pointer_1;
char* pointer_2;
int primer_1;
int primer_2;
HP_DynamicStringArray* barcodes = new HP_DynamicStringArray();
HPDynamicIntArray* up_down = new HPDynamicIntArray();
char* temp = new char[21];
for(int x = 0; x < string_array->GetNumberOfStrings(); x++)
	{
	memset(temp, 0, 21);
	int index = 0;
	while(index < 4 && !strstr(string_array->GetStringPointer(x), primers[index]))
		index++;
	if(index < 4) // we  found a match
		{
		pointer_1 = strstr(string_array->GetStringPointer(x), primers[index]);
		// See if we can find the match for the primer pair
		pointer_2 = strstr(string_array->GetStringPointer(x), primers[index+4]);
		primer_1 = pointer_1 - string_array->GetStringPointer(x);
		primer_2 = pointer_2 - string_array->GetStringPointer(x);
		if(pointer_2 == 0) // could not find a match
			{
				one_primer++;
			}
		else
			{
			if(primer_2 > primer_1) // the order of primers is correct
				{
				two_primers++;
				if((primer_2 - primer_1 - primer_length) == insert_length)
					{
						correct_length++;
						strncpy(temp, (pointer_1+primer_length), insert_length);
						if(index < 2)
							up_down->AddInt(1);
						else
							up_down->AddInt(0);
						if((index == 1) || (index == 3))
							barcodes->AddString(readfile->ReverseComplement(temp));
						else
							barcodes->AddString(temp);


					}
				else
					incorrect_length++;
				}
			}

		}
	else
		{
		no_primers++;
		}
	}

cout << "two_primers: " << two_primers << endl;
cout << "correct_length: " << correct_length << endl;
cout << "one_primer: " << one_primer << endl;
cout << "incorrect_length: " << incorrect_length << endl;
cout << "no_primers: " << no_primers << endl;

// Save the strings from the array

char* outfilename = new char[1024];

//	Change the string below to a filename and path of your choice

fstream outfile;
outfile.open(outfile_filepath, fstream::out | fstream::binary);
char* temp_2 = new char[30];
for(int x = 0; x < barcodes->GetNumberOfStrings(); x++)
{
	memset(temp_2, 0, 30);
	sprintf(temp_2, "%s%c%d", barcodes->GetStringPointer(x), '\t', up_down->GetEntry(x));
	outfile << temp_2 << endl;
}
outfile.close();

cout << "done" << endl;

return 0;
}
