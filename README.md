 const1.h - A collection of very handy tools for SIMD vector development which
            most of my projects on this repo rely on

 Note: you must include simd_c.c into project for ISINF, ISNAN, 
       ISNORM, FIXNORM to work properly

SIMD's simd::isinf simd::isnan is stil broken even in the latest Xcode 
(16.2), so I built my own from the ground up, would be nice if Apple would
fix this on compiler level, but for now here is a replacement, along
with other many handy tools I accumulated devleoping projects, including 
ZArray, which I am looking for support on compiler level as well, able 
to quickly reference an array of float vectors by an integer vector 
returning a float vector with values corresponding to the indices
of the index vector.

This bug in SIMD isinf, isnan has been persistent through many of your 
XCode releases already, and management seem to out to lunch about it, 
so get to work, Apple and stop slacking!

Copyright 2022-2025 Dmitry Boldyrev
All rights reserved.

Report bugs, issues to my email: subband@gmail.com  
