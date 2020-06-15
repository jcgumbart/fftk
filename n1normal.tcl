        # handle the n=1 case 'normally'

        # find some information about B atom
        set bInd $bondlistA
        set bSel [atomselect top "index $bInd"]
        set bondlistB [lindex [$bSel getbonds] 0]

        # find a valid C atom
        set cInd {}
        foreach ele $bondlistB {
            # make sure that C != A, and is non-linear
            if { $ele == $aInd } {
                continue
            } elseif { [expr {abs([measure angle [list $aInd $bInd $ele]])}] >= 178.0 } {
                continue
            } else {
                set cInd $ele
            }
        }

        # make an atom selection of C atom
        set cSel [atomselect top "index $cInd"]

        # get gnames
        set aGname [lindex $gnames $aInd]
        set bGname [lindex $gnames $bInd]
        set cGname [lindex $gnames $cInd]

        if { $class eq "donor" } {
            # donor
            puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" x $aGname 1.0 $bGname 90.00 $cGname dih]
            puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" Ow $aGname rAH x 90.00 $bGname 180.00]
            puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" H1w Ow 0.9572 $aGname 127.74 x 0.00]
            puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s\n" H2w Ow 0.9572 $aGname 127.74 x 180.00]
            puts $outfile "rAH  2.0"
            puts $outfile "dih  0.0"
        } else {
            # acceptor
            # call helper function to find probe position
            set probePos [::ForceFieldToolKit::Gaussian::placeProbe $aSel]
            # note that Gaussian doesn't like 180 degree angles, so we have to be a little clever
            # make some measurements of probe position
            set mAng [::QMtool::bond_angle $probePos [measure center $aSel] [measure center $cSel]]
            set mDih [::QMtool::dihed_angle $probePos [measure center $aSel] [measure center $cSel] [measure center $bSel]]

            puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" H1w $aGname rAH $cGname [format %3.2f $mAng] $bGname [format %3.2f $mDih]]
            puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" x H1w 1.0 $aGname 90.00 $cGname 0.00]
            puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" Ow H1w 0.9527 x 90.00 $aGname 180.00]
            puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s\n" H2w Ow 0.9527 H1w 104.52 x dih]
            puts $outfile "rAH  2.0"
            puts $outfile "dih  0.0"
        }

        # clean up and return
        $aSel delete; $bSel delete; $cSel delete
        return			
