		#=======#
		# N > 1 #
		#=======#

		# find some information about atom B
		set bInd [lindex $bondlistA 0]
		set bSel [atomselect top "index $bInd"]
		set cInd [lindex $bondlistA 1]
		set cSel [atomselect top "index $cInd"]

		# find gnames
		set aGname [lindex $gnames $aInd]
		set bGname [lindex $gnames $bInd]
		set cGname [lindex $gnames $cInd]

		# test if C is valid choice
		set abcAng [expr {abs([measure angle [list $aInd $bInd $cInd]])}]
		if { $abcAng < 2 || $abcAng > 178 } { set validC 0 } else {	set validC 1 }
		unset abcAng

		# find probe coords
		set probePos [::ForceFieldToolKit::Gaussian::placeProbe $aSel]
		set mAng [::QMtool::bond_angle $probePos [measure center $aSel] [measure center $bSel]]

		if { !$validC } {
			# if C is invalid, ABC are linear and we need a second dummy atom in lieu of the original C atom
			set cGname "x2"
			set x2Pos [coordtrans [trans center [measure center $aSel] axis [vecsub [measure center $cSel] [measure center $bSel]] 180.0 deg] $probePos]
			set mDih [::QMtool::dihed_angle $probePos [measure center $aSel] [measure center $bSel] $x2Pos]
			puts $outfile [format "%3s  %.4f  %.4f  %.4f" x2 [lindex $x2Pos 0] [lindex $x2Pos 1] [lindex $x2Pos 2]]
		} else {
			# C is valid, we can use it to define the dihedral
			set mDih [::QMtool::dihed_angle $probePos [measure center $aSel] [measure center $bSel] [measure center $cSel]]	
		}
		
		if { $class eq "donor" } {
			# donor
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" Ow $aGname rAH $bGname [format %3.2f $mAng] $cGname [format %3.2f $mDih]]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" x Ow 1.0 $aGname 90.00 $bGname dih]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" H1w Ow 0.9572 $aGname 127.74 x 0.00]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s\n" H2w Ow 0.9572 $aGname 127.74 x 180.00]
            puts $outfile "rAH  2.0"
            puts $outfile "dih  0.0"

		} else {
			# acceptor
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" H1w $aGname rAH $bGname [format %3.2f $mAng] $cGname [format %3.2f $mDih]]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" x H1w 1.0 $aGname 90.00 $bGname 0.00]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" Ow H1w 0.9572 x 90.00 $aGname 180.00]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s\n" H2w Ow 0.9572 H1w 104.52 x dih]
			puts $outfile "rAH  2.0"
			puts $outfile "dih  0.0"
		}

		# clean up atomselections and return
		$aSel delete; $bSel delete; $cSel delete
		return
