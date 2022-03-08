/** 
* Copyright (c) 2013, Laboratory for Biocomputing and Informatics, Boston University
*  All rights reserved.
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions are met:
*   + Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
*   + Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*   + This source code may not be used in any program that is sold, any
*    derivative work must allow free distribution and modification.
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
*  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED
*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
*  DISCLAIMED. IN NO EVENT SHALL THE LABORATORY FOR BIOCOMPUTING AND INFORMATICS 
*  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
OR CONSEQUENTIAL 
*  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
*  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
TORT
*  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
*  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
        
//Uses 75 operations, counting |, &, ^, and +

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define wordsize 64

int bitwise_alignment_m0_x1_g1(char * s1, char* s2, int words) {
	
	int debug = 1;
	int wordsizem1 = wordsize - 1;
	int wordsizem2 = wordsize - 2;
	int i, j, k, n, m;
	unsigned long long int bitmask;
	unsigned long long int Matches;
	unsigned long long int **matchvec;
	unsigned long long int *matchvecmem;
	unsigned long long int *matchv;
	unsigned long long int NotMatches;
	unsigned long long int all_ones = ~0x0000000000000000;
	unsigned long long int one = 0x0000000000000001;
	unsigned long long int sixtythreebitmask = ~((unsigned long long int) 1 << wordsizem1);
	
	unsigned long long int DVpos1shift, DVzeroshift, DVneg1shift;	unsigned long long int INITpos1s, INITzeros;
	unsigned long long int DHpos1, DHzero, DHneg1;
	unsigned long long int INITpos1sprevbit, INITzerosprevbit;
	unsigned long long int OverFlow0, OverFlow1;
	
	unsigned long long int RemainDHneg1;
	unsigned long long int DHneg1tozero;
	unsigned long long int DVpos1shiftorMatch;
	unsigned long long int DVnot1to1shiftorMatch;
	unsigned long long int DHpos1orMatch;
	
	//unsigned long long int DHneg3toDHneg2;
	//unsigned long long int DHneg1toDHzer0;
	//unsigned long long int DHpos1toDHpos2;
	//unsigned long long int DHpos3toDHpos4;
	//unsigned long long int DHpos5toDHpos6;
	//unsigned long long int DHneg1topos2;
	//unsigned long long int DHpos3topos6;
	unsigned long long int notDHbits4;
	unsigned long long int DHbits4notDHbits2;
	unsigned long long int DHbits4DHbits2;
	unsigned long long int notDHbits2;
	unsigned long long int notDHbits4DHbits2;
	unsigned long long int notDHbits4notDHbits2;
	unsigned long long int notDHbits1;
	unsigned long long int DVbits4DVbits2;
	unsigned long long int notDVbits4;
	unsigned long long int notDVbits4DVbits2;
	unsigned long long int notDVbits2;
	unsigned long long int notDVbits1;
	unsigned long long int DVbits4notDVbits2;
	unsigned long long int notDVbits4notDVbits2;
	unsigned long long int * DVDHbits4;
	unsigned long long int * DVDHbits2;
	unsigned long long int * DVDHbits1;
	unsigned long long int * DVDHbits8;
	unsigned long long int carry0;
	unsigned long long int carry1;
	unsigned long long int carry2;
	unsigned long long int DVDHbithighcomp;
	unsigned long long int compDHneg1tozero;
	unsigned long long int DVbits1;
	unsigned long long int DVbits4;
	unsigned long long int DVbits2;
	unsigned long long int DHbits2;
	unsigned long long int DHbits1;
	unsigned long long int DHbits4;
	unsigned long long int DHpos1orDHzero;
	unsigned long long int tempbit4,tempbit2,tempbit1;
	unsigned long long int sumbit1,sumbit2,sumbit4;
	unsigned long long int bit4,bit2,bit1;
	
	unsigned long long int sum;
	unsigned long long int highone = one << wordsizem1;
	unsigned long long int nexthighone = one << wordsizem2;
	unsigned long long int initval;
	char * iterate;
	int counter = 0, w = 0;
	int score;
	
	n = strlen (s1);
	m = strlen (s2);
	
	//*************************encode match strings A C G T N for string1
	//loop through string1 and store bits in matchA, matchC, etc.
	//position zero corresponds to column one in the score matrix, i.e., first character
	//so we start with i = 0 and bitmask = 1
	
	bitmask = one;
	matchvec = (unsigned long long int **) calloc(256, sizeof(unsigned long long int *));
	matchvecmem = (unsigned long long int *) calloc(words * 256, sizeof(unsigned long long int));
	for (i = 0 ; i < 256; ++i)
	matchvec[i] = &matchvecmem[i*words];
	for (iterate = s1, i = 0; i < n; ++i, ++iterate)
	{
		matchvec[(*iterate)][w] |= bitmask;
		bitmask <<= 1; ++counter;
		if (counter == 63)
		{
			counter = 0;
			w++;
			bitmask = one;
		}
	}
	DVDHbits1 = (unsigned long long int * )calloc(words, sizeof(unsigned long long int));
	DVDHbits2 = (unsigned long long int * )calloc(words, sizeof(unsigned long long int));
	DVDHbits4 = (unsigned long long int * )calloc(words, sizeof(unsigned long long int));
	DVDHbits8 = (unsigned long long int * )calloc(words, sizeof(unsigned long long int));
	
	
	
	//intialize top row (penalty for initial gap)
	DHneg1 = all_ones;
	DHzero = DHpos1 = 0;
	DHpos1orDHzero = (DHpos1|DHzero);
	for (i = 0; i < words; ++i){
		DVDHbits1[i] = DHzero;
		DVDHbits2[i] = DHpos1orDHzero;
		DVDHbits4[i] = DHpos1orDHzero;
		DVDHbits8[i] = 0;
	}
	
	//recursion
	for (i = 0, iterate = s2; i < m; ++i, ++iterate)
	{
		bit1 = bit2 = bit4 = 0;
		//initialize OverFlow
		OverFlow0 = OverFlow1 = 0;
		INITpos1sprevbit = INITzerosprevbit = 0;
		
		matchv = matchvec[*iterate];
		for (j = 0; j < words; ++j){
			DHbits4 = DVDHbits4[j];
			DHbits2 = DVDHbits2[j];
			DHbits1 = DVDHbits1[j];
			
			Matches = *matchv;
			matchv++;
			//Complement Matches
			NotMatches = ~Matches;
			notDHbits4 = (~DHbits4);
			notDHbits2 = ~DHbits2;
			DHbits4DHbits2 = DHbits4 & DHbits2;
			notDHbits4notDHbits2 = notDHbits4 & notDHbits2;
			notDHbits1 = ~DHbits1;
			DHneg1 = notDHbits4notDHbits2 & notDHbits1 & sixtythreebitmask;
			DHzero = DHbits4DHbits2 & DHbits1;
			//Finding the vertical values
			//Find 1s
			INITpos1s = DHneg1 & Matches;
			sum = (INITpos1s + DHneg1) + OverFlow0;
			DVpos1shift = ((sum  ^ DHneg1) ^ INITpos1s) & sixtythreebitmask;
			OverFlow0 = sum >> wordsizem1;
			
			//set RemainingDHneg1
			RemainDHneg1 = DHneg1 ^ (INITpos1s);
			//combine 1s and Matches
			DVpos1shiftorMatch = DVpos1shift | Matches;
			
			//set DVnot1to1shiftorMatch
			DVnot1to1shiftorMatch = ~(DVpos1shiftorMatch);
			DVbits1 = DVnot1to1shiftorMatch;
			DVbits2 = DVpos1shiftorMatch;
			DVbits4 = 0;
			carry0 = (DHbits1&DVbits1);
			carry1 = ((DHbits2&DVbits2) | ((DHbits2^DVbits2)&carry0));
			carry2 = ((DHbits4&DVbits4) | ((DHbits4^DVbits4)&carry1));
			sumbit1 = (DHbits1^DVbits1);
			sumbit2 = ((DHbits2^DVbits2)^carry0);
			sumbit4 = ((DHbits4^DVbits4)^carry1);
			DVDHbithighcomp = ~sumbit4;
			sumbit1 &= DVDHbithighcomp;
			sumbit2 &= DVDHbithighcomp;
			sumbit4 &= DVDHbithighcomp;
			tempbit4 = (sumbit4 & nexthighone) >> wordsizem2;
			tempbit2 = (sumbit2 & nexthighone) >> wordsizem2;
			tempbit1 = (sumbit1 & nexthighone) >> wordsizem2;
			sumbit1 <<= 1;
			sumbit2 <<= 1;
			sumbit4 <<= 1;
			sumbit4 |= bit4;
			sumbit2 |= bit2;
			sumbit1 |= bit1;
			bit4 = tempbit4;
			bit2 = tempbit2;
			bit1 = tempbit1;
			DVbits1 = sumbit1;
			DVbits2 = sumbit2;
			DVbits4 = sumbit4;
			
			//Finding the horizontal values
			//Remove matches from DH values except 1
			//combine 1s and Matches
			DHpos1orMatch = DHpos1| Matches;
			//group -1topos0
			DHneg1tozero = (DHneg1 | DHzero) & NotMatches;
			compDHneg1tozero = all_ones ^ DHneg1tozero;
			//Representing the range -1 to 0 as 0;
			DHbits1 |= DHneg1tozero;
			DHbits2 |= DHneg1tozero;
			DHbits4 |= DHneg1tozero;
			//Representing Matches as the packed value of 1 which is -2
			DHbits1 &= NotMatches;
			DHbits2 |= Matches;
			DHbits4 |= Matches;
			carry0 = (DHbits1&sumbit1);
			carry1 = ((DHbits2&sumbit2) | ((DHbits2^sumbit2)&carry0));
			carry2 = ((DHbits4&sumbit4) | ((DHbits4^sumbit4)&carry1));
			sumbit1 = (DHbits1^sumbit1);
			sumbit2 = ((DHbits2^sumbit2)^carry0);
			sumbit4 = ((DHbits4^sumbit4)^carry1);
			//Convert everything that is positive to 0
			sumbit1 &= sumbit4;
			sumbit2 &= sumbit4;
			sumbit4 &= sumbit4;
			DVDHbits1[j] = sumbit1;
			DVDHbits2[j] = sumbit2;
			DVDHbits4[j] = sumbit4;
			
			
		}
	}
	//find scores in last row
	
	score = -1 * m;
	
	bitmask = one;
	for (j = 0; j < words; ++j){
		for (i = j*wordsizem1; i < (j + 1)*wordsizem1 && i < n; ++i)
		{
			score -= (DVDHbits1[j] & bitmask) * 1 + (DVDHbits2[j] & bitmask) * 2 + -(DVDHbits4[j] & bitmask) * 4 + 1;
			DVDHbits1[j] >>= 1;
			DVDHbits2[j] >>= 1;
			DVDHbits4[j] >>= 1;
			
			
		}
		
	}
	free(DVDHbits1);
	free(DVDHbits2);
	free(DVDHbits4);
	free(DVDHbits8);
	
	free(matchvecmem);
	free(matchvec);
	
	return score;
}

	
