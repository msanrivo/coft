

    LEAPSECONDS KERNEL VERSION NAIF0010

The file naif0010.tls is a unix-style text file. It is suitable for use
on all unix boxes, including PC/linux and MAC/OSX machines.

For PCs running Windows, use naif0010.tls.pc.

Use of one of these files is required for all SPICE computations
dealing with times on or after July 01, 2012 00:00:00 UTC if you want
accurate conversions between UTC and Ephemeris Time (a.k.a. TDB).
Failure to use naif0010 for times after this epoch will result
in a one second time conversion error.

You may begin using naif0010 in place of naif0009 right now; there
is no need to wait until July 01, 2012. Time conversions for epochs
prior to July 01, 2012 will be the same using either naif0009 or
naif0010.


