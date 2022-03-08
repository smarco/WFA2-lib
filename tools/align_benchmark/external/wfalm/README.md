# wfalm
Refinements of the WFA alignment algorithm with better complexity

## Introduction

This repository contains implementations of the [WFA algorithm](https://academic.oup.com/bioinformatics/article/37/4/456/5904262) as well as two variations on it. The first reduces the memory complexity from O(s^2) to O(s^3/2). The second reduces the time complexity from O(sN) to O(s^2 + N) using suffix trees. The former is eminently practical, but the latter has large constants that make it take more time in the majority of common use cases.

## Installation

Because many users will not need the suffix tree algorithm, the library is implemented in two headers. The non-suffix tree algorithm is implemented in `wfa_lm.hpp`. It is single-header and has no external dependencies beyond the 2011 C++ standard library. It can simply be included into a new project.

The suffix tree algorithm is implemented in `wfa_lm_st.hpp`. It is also single-header, but it depends on [SDSL v.3](https://github.com/xxsds/sdsl-lite) for its implementation of a suffix tree. This dependency can be pulled in with this command:

	git submodule update --init --recursive
	
After that, the headers in `external/sdsl-lite/include` will also need to be included in the new project.