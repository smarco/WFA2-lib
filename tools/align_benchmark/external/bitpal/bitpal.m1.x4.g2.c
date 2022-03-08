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
        
//Uses 154 operations, counting |, &, ^, and +

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define wordsize 64

int bitwise_alignment_m1_x4_g2(char * s1, char* s2, int words) {
	
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
	
	unsigned long long int DVpos3shift, DVpos2shift, DVpos1shift, DVzeroshift, DVneg1shift, DVneg2shift;	unsigned long long int INITpos3s, INITpos2s, INITpos1s, INITzeros, INITneg1s;
	unsigned long long int DHpos3, DHpos2, DHpos1, DHzero, DHneg1, DHneg2;
	unsigned long long int INITpos3sprevbit, INITpos2sprevbit, INITpos1sprevbit, INITzerosprevbit, INITneg1sprevbit;
	unsigned long long int OverFlow0, OverFlow1, OverFlow2, OverFlow3, OverFlow4;
	
	unsigned long long int DVpos2shiftNotMatch, DVpos1shiftNotMatch, DVzeroshiftNotMatch, DVneg1shiftNotMatch;
	unsigned long long int RemainDHneg2;
	unsigned long long int DHneg2toneg2;
	unsigned long long int DVpos3shiftorMatch;
	unsigned long long int DVnot3to4shiftorMatch;
	unsigned long long int DHpos3orMatch;
	
	//unsigned long long int DHneg3toDHneg2;
	//unsigned long long int DHneg1toDHzer0;
	//unsigned long long int DHpos1toDHpos2;
	//unsigned long long int DHpos3toDHpos4;
	//unsigned long long int DHpos5toDHpos6;
	//unsigned long long int DHneg1topos2;
	//unsigned long long int DHpos3topos6;
	unsigned long long int notDHbits4;
	unsigned long long int notDHbits8DHbits4;
	unsigned long long int notDHbits8;
	unsigned long long int DHbits8DHbits4DHbits2;
	unsigned long long int DHbits8notDHbits4notDHbits2;
	unsigned long long int notDHbits8DHbits4DHbits2;
	unsigned long long int DHbits8notDHbits4DHbits2;
	unsigned long long int notDHbits8notDHbits4DHbits2;
	unsigned long long int notDHbits1;
	unsigned long long int notDHbits8notDHbits4notDHbits2;
	unsigned long long int notDHbits2;
	unsigned long long int DHbits8notDHbits4;
	unsigned long long int notDHbits8notDHbits4;
	unsigned long long int notDHbits8DHbits4notDHbits2;
	unsigned long long int DHbits8DHbits4;
	unsigned long long int DHbits8DHbits4notDHbits2;
	unsigned long long int DVbits8DVbits4notDVbits2;
	unsigned long long int notDVbits8DVbits4notDVbits2;
	unsigned long long int notDVbits8DVbits4DVbits2;
	unsigned long long int notDVbits4;
	unsigned long long int DVbits8DVbits4;
	unsigned long long int notDVbits8DVbits4;
	unsigned long long int DVbits8notDVbits4;
	unsigned long long int notDVbits8notDVbits4DVbits2;
	unsigned long long int notDVbits2;
	unsigned long long int notDVbits1;
	unsigned long long int notDVbits8;
	unsigned long long int notDVbits8notDVbits4notDVbits2;
	unsigned long long int notDVbits8notDVbits4;
	unsigned long long int DVbits8notDVbits4notDVbits2;
	unsigned long long int DVbits8DVbits4DVbits2;
	unsigned long long int DVbits8notDVbits4DVbits2;
	unsigned long long int * DVDHbits8;
	unsigned long long int * DVDHbits4;
	unsigned long long int * DVDHbits2;
	unsigned long long int * DVDHbits1;
	unsigned long long int * DVDHbits16;
	unsigned long long int carry0;
	unsigned long long int carry1;
	unsigned long long int carry2;
	unsigned long long int carry3;
	unsigned long long int DVDHbithighcomp;
	unsigned long long int compDHneg2toneg2;
	unsigned long long int DVpos2shiftorDVpos3shiftorMatch;
	unsigned long long int DVzeroshiftorDVpos1shift;
	unsigned long long int DVbits4;
	unsigned long long int DVbits1;
	unsigned long long int DVbits8;
	unsigned long long int DVbits2;
	unsigned long long int DHbits2;
	unsigned long long int DHpos2orDHpos1;
	unsigned long long int DHbits8;
	unsigned long long int DHbits1;
	unsigned long long int DHpos2orDHpos1orDHzeroorDHneg1;
	unsigned long long int DHzeroorDHneg1;
	unsigned long long int DHpos3orDHpos2orDHpos1orDHzeroorDHneg1;
	unsigned long long int DHbits4;
	unsigned long long int tempbit8,tempbit4,tempbit2,tempbit1;
	unsigned long long int sumbit1,sumbit2,sumbit4,sumbit8;
	unsigned long long int bit8,bit4,bit2,bit1;
	
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
	DVDHbits16 = (unsigned long long int * )calloc(words, sizeof(unsigned long long int));
	
	
	
	//intialize top row (penalty for initial gap)
	DHneg2 = all_ones;
	DHneg1 = DHzero = DHpos1 = DHpos2 = DHpos3 = 0;
	DHzeroorDHneg1 = (DHzero|DHneg1);
	DHpos2orDHpos1 = (DHpos2|DHpos1);
	DHpos2orDHpos1orDHzeroorDHneg1 = (DHpos2orDHpos1|DHzeroorDHneg1);
	DHpos3orDHpos2orDHpos1orDHzeroorDHneg1 = (DHpos3|DHpos2orDHpos1orDHzeroorDHneg1);
	for (i = 0; i < words; ++i){
		DVDHbits1[i] = DHneg1 | DHpos1 | DHpos3;
		DVDHbits2[i] = DHzeroorDHneg1 | DHpos3;
		DVDHbits4[i] = DHpos2orDHpos1orDHzeroorDHneg1;
		DVDHbits8[i] = DHpos3orDHpos2orDHpos1orDHzeroorDHneg1;
		DVDHbits16[i] = 0;
	}
	
	//recursion
	for (i = 0, iterate = s2; i < m; ++i, ++iterate)
	{
		bit1 = bit2 = bit4 = bit8 = 0;
		//initialize OverFlow
		OverFlow0 = OverFlow1 = OverFlow2 = OverFlow3 = OverFlow4 = 0;
		INITpos3sprevbit = INITpos2sprevbit = INITpos1sprevbit = INITzerosprevbit = INITneg1sprevbit = 0;
		
		matchv = matchvec[*iterate];
		for (j = 0; j < words; ++j){
			DHbits8 = DVDHbits8[j];
			DHbits4 = DVDHbits4[j];
			DHbits2 = DVDHbits2[j];
			DHbits1 = DVDHbits1[j];
			
			Matches = *matchv;
			matchv++;
			//Complement Matches
			NotMatches = ~Matches;
			notDHbits8 = (~DHbits8);
			notDHbits4 = ~DHbits4;
			notDHbits8notDHbits4 = notDHbits8 & notDHbits4;
			notDHbits2 = ~DHbits2;
			notDHbits8notDHbits4notDHbits2 = notDHbits8notDHbits4 & notDHbits2;
			notDHbits1 = ~DHbits1;
			DHneg2 = notDHbits8notDHbits4notDHbits2 & notDHbits1 & sixtythreebitmask;
			//Finding the vertical values
			//Find 3s
			INITpos3s = DHneg2 & Matches;
			sum = (INITpos3s + DHneg2) + OverFlow0;
			DVpos3shift = ((sum  ^ DHneg2) ^ INITpos3s) & sixtythreebitmask;
			OverFlow0 = sum >> wordsizem1;
			
			//set RemainingDHneg2
			RemainDHneg2 = DHneg2 ^ (INITpos3s);
			//combine 3s and Matches
			DVpos3shiftorMatch = DVpos3shift | Matches;
			
			//Find 2s
			INITpos2s = (DHneg1 & DVpos3shiftorMatch) ;
			initval = ((INITpos2s << 1) | INITpos2sprevbit);
			INITpos2sprevbit = (initval >> wordsizem1);
			initval &= sixtythreebitmask;
			sum = initval + RemainDHneg2 + OverFlow1;
			DVpos2shift = (sum  ^ RemainDHneg2) & NotMatches;
			OverFlow1= sum >> wordsizem1;
			//Find 1s
			INITpos1s = (DHzero & DVpos3shiftorMatch) | (DHneg1 & DVpos2shift);
			initval = ((INITpos1s << 1) | INITpos1sprevbit);
			INITpos1sprevbit = (initval >> wordsizem1);
			initval &= sixtythreebitmask;
			sum = initval + RemainDHneg2 + OverFlow2;
			DVpos1shift = (sum  ^ RemainDHneg2) & NotMatches;
			OverFlow2= sum >> wordsizem1;
			//Find 0s
			INITzeros = (DHpos1 & DVpos3shiftorMatch) | (DHzero & DVpos2shift)| (DHneg1 & DVpos1shift);
			initval = ((INITzeros << 1) | INITzerosprevbit);
			INITzerosprevbit = (initval >> wordsizem1);
			initval &= sixtythreebitmask;
			sum = initval + RemainDHneg2 + OverFlow3;
			DVzeroshift = (sum  ^ RemainDHneg2) & NotMatches;
			OverFlow3= sum >> wordsizem1;
			//Find -1s
			INITneg1s = (DHpos2 & DVpos3shiftorMatch) | (DHpos1 & DVpos2shift)| (DHzero & DVpos1shift)| (DHneg1 & DVzeroshift);
			initval = ((INITneg1s << 1) | INITneg1sprevbit);
			INITneg1sprevbit = (initval >> wordsizem1);
			initval &= sixtythreebitmask;
			sum = initval + RemainDHneg2 + OverFlow4;
			DVneg1shift = (sum  ^ RemainDHneg2) & NotMatches;
			OverFlow4= sum >> wordsizem1;
			//set DVnot3to4shiftorMatch
			DVnot3to4shiftorMatch = ~(DVpos3shiftorMatch | DVpos2shift | DVpos1shift | DVzeroshift | DVneg1shift);
			DVpos2shiftorDVpos3shiftorMatch = (DVpos2shift|DVpos3shiftorMatch);
			DVzeroshiftorDVpos1shift = (DVzeroshift|DVpos1shift);
			DVbits1 = DVpos3shiftorMatch | DVpos1shift | DVneg1shift;
			DVbits2 = DVzeroshiftorDVpos1shift;
			DVbits4 = DVpos2shiftorDVpos3shiftorMatch;
			DVbits8 = 0;
			carry0 = (DHbits1&DVbits1);
			carry1 = ((DHbits2&DVbits2) | ((DHbits2^DVbits2)&carry0));
			carry2 = ((DHbits4&DVbits4) | ((DHbits4^DVbits4)&carry1));
			carry3 = ((DHbits8&DVbits8) | ((DHbits8^DVbits8)&carry2));
			sumbit1 = (DHbits1^DVbits1);
			sumbit2 = ((DHbits2^DVbits2)^carry0);
			sumbit4 = ((DHbits4^DVbits4)^carry1);
			sumbit8 = ((DHbits8^DVbits8)^carry2);
			DVDHbithighcomp = ~sumbit8;
			sumbit1 &= DVDHbithighcomp;
			sumbit2 &= DVDHbithighcomp;
			sumbit4 &= DVDHbithighcomp;
			sumbit8 &= DVDHbithighcomp;
			tempbit8 = (sumbit8 & nexthighone) >> wordsizem2;
			tempbit4 = (sumbit4 & nexthighone) >> wordsizem2;
			tempbit2 = (sumbit2 & nexthighone) >> wordsizem2;
			tempbit1 = (sumbit1 & nexthighone) >> wordsizem2;
			sumbit1 <<= 1;
			sumbit2 <<= 1;
			sumbit4 <<= 1;
			sumbit8 <<= 1;
			sumbit8 |= bit8;
			sumbit4 |= bit4;
			sumbit2 |= bit2;
			sumbit1 |= bit1;
			bit8 = tempbit8;
			bit4 = tempbit4;
			bit2 = tempbit2;
			bit1 = tempbit1;
			DVbits1 = sumbit1;
			DVbits2 = sumbit2;
			DVbits4 = sumbit4;
			DVbits8 = sumbit8;
			
			//Finding the horizontal values
			//Remove matches from DH values except 3
			//combine 3s and Matches
			DHpos3orMatch = DHpos3| Matches;
			//group -2topos1
			DHneg2toneg2 = (DHneg1) & NotMatches;
			compDHneg2toneg2 = all_ones ^ DHneg2toneg2;
			//Representing the range -2 to -2 as -2;
			DHbits1 &= compDHneg2toneg2;
			DHbits2 &= compDHneg2toneg2;
			DHbits4 &= compDHneg2toneg2;
			DHbits8 &= compDHneg2toneg2;
			//Representing Matches as the packed value of 3 which is -5
			DHbits1 |= Matches;
			DHbits2 |= Matches;
			DHbits4 &= NotMatches;
			DHbits8 |= Matches;
			carry0 = (DHbits1&sumbit1);
			carry1 = ((DHbits2&sumbit2) | ((DHbits2^sumbit2)&carry0));
			carry2 = ((DHbits4&sumbit4) | ((DHbits4^sumbit4)&carry1));
			carry3 = ((DHbits8&sumbit8) | ((DHbits8^sumbit8)&carry2));
			sumbit1 = (DHbits1^sumbit1);
			sumbit2 = ((DHbits2^sumbit2)^carry0);
			sumbit4 = ((DHbits4^sumbit4)^carry1);
			sumbit8 = ((DHbits8^sumbit8)^carry2);
			//Convert everything that is positive to 0
			sumbit1 &= sumbit8;
			sumbit2 &= sumbit8;
			sumbit4 &= sumbit8;
			sumbit8 &= sumbit8;
			DVDHbits1[j] = sumbit1;
			DVDHbits2[j] = sumbit2;
			DVDHbits4[j] = sumbit4;
			DVDHbits8[j] = sumbit8;
			
			
		}
	}
	//find scores in last row
	
	score = -2 * m;
	
	bitmask = one;
	for (j = 0; j < words; ++j){
		for (i = j*wordsizem1; i < (j + 1)*wordsizem1 && i < n; ++i)
		{
			score -= (DVDHbits1[j] & bitmask) * 1 + (DVDHbits2[j] & bitmask) * 2 + (DVDHbits4[j] & bitmask) * 4 + -(DVDHbits8[j] & bitmask) * 8 + 2;
			DVDHbits1[j] >>= 1;
			DVDHbits2[j] >>= 1;
			DVDHbits4[j] >>= 1;
			DVDHbits8[j] >>= 1;
			
			
		}
		
	}
	free(DVDHbits1);
	free(DVDHbits2);
	free(DVDHbits4);
	free(DVDHbits8);
	free(DVDHbits16);
	
	free(matchvecmem);
	free(matchvec);
	
	return score;
}

