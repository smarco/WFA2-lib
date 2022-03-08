This repo implements the [Landau-Vishkin algorithm][lv89] to compute the edit
distance between two strings. This is a fast method for highly similar strings.
[The actual implementation](lv89.c) follows a simplified [WFA][WFA] formulation
rather than the original formulation. It also learns a performance trick from
WFA. For a pair of ~5Mb HLA sequences with ~123k edits, the implementation here
can find the result in 71 seconds, faster than [edlib][edlib]. Edlib will be
faster for more divergent sequences.

[lv89]: https://doi.org/10.1016/0196-6774(89)90010-2
[edlib]: https://github.com/Martinsos/edlib
[WFA]: https://github.com/smarco/WFA
