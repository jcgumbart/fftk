#===========================================================================================================
#===========================================================================================================
# ffTK QM Module Support for Psi4
#===========================================================================================================
#===========================================================================================================
namespace eval ::ForceFieldToolKit::Psi4 {

}
#===========================================================================================================
# GEOMETRY OPTIMIZATION
#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::resetDefaultsGeomOpt {} {
    # resets the Psi4 settings to default

    # reset to default
    set ::ForceFieldToolKit::GeomOpt::qmProc 1
    set ::ForceFieldToolKit::GeomOpt::qmMem 1
    set ::ForceFieldToolKit::GeomOpt::qmCharge 0
    set ::ForceFieldToolKit::GeomOpt::qmMult 1
    set ::ForceFieldToolKit::GeomOpt::qmRoute "mp2/6-31G*"

}
#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::writeComGeomOpt { molID com qmProc qmMem qmCharge qmMult qmRoute } {
    # write input file for geometry optimization

    # assign atom names and gather x,y,z for output com file
    set Gnames {}
    set atom_info {}

    for {set i 0} {$i < [molinfo $molID get numatoms]} {incr i} {
        set temp [atomselect $molID "index $i"]
        lappend atom_info [list [$temp get element] [$temp get x] [$temp get y] [$temp get z]]
        lappend Gnames [$temp get element][expr $i+1]
        $temp delete
    }

    # open the output com file
    set outfile [open $com w]

    # start python code
    puts $outfile "import psi4"
    puts $outfile ""
    puts $outfile "psi4.set_output_file(\"output.out\", False)"
    puts $outfile "psi4.set_memory(\"$qmMem GB\")"
    puts $outfile ""
    puts $outfile {molecule = psi4.geometry("""}

    # write the coordinates
    foreach atom_entry $atom_info {
        puts $outfile "[lindex $atom_entry 0] [format %16.8f [lindex $atom_entry 1]] [format %16.8f [lindex $atom_entry 2]] [format %16.8f [lindex $atom_entry 3]]"
   }

    puts $outfile {""")}
    puts $outfile ""
    puts $outfile "molecule.set_multiplicity($qmMult)"
    puts $outfile "molecule.set_molecular_charge($qmCharge)"
    puts $outfile "psi4.optimize(\"$qmRoute\", molecule=molecule)"

    # empty line to terminate
    puts $outfile ""

    # clean up
   close $outfile

}

#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::readOutGeomOpt { molId logFile } {
    # load output file from geometry optimization

    # load the coordinates from Psi4 log
    set inFile [open $logFile r]

    while { ![eof $inFile] } {
        set inLine [string trim [gets $inFile]]
        if { $inLine eq "OPTKING Finished Execution" } {
            # burn the coord header
            for {set i 0} {$i < 8} {incr i} { gets $inFile }
            # read coordinates
            set coords {}
            while { ![regexp {^-*$} [set inLine [string trim [gets $inFile]]]] } {
                lappend coords [lrange $inLine 1 4]
            }
            set coords [::ForceFieldToolKit::SharedFcns::LonePair::addLPCoordinate $coords]

            # add a new frame, set the coords
            animate dup $molId
            set temp [atomselect $molId all]
            $temp set {x y z} $coords
            $temp delete
            unset coords
        } else {
            continue
        }
    }

    # clean up
    close $inFile

}
#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::writePDBGeomOpt { psf pdb logFile optPdb } {
    # write the optimized structure to as new PDB file

    # load the pdb and load the coords from the log file
    set molId [::ForceFieldToolKit::SharedFcns::LonePair::initFromPSF $psf]
    mol addfile $pdb

    # NEW 02/01/2019: Checking "element" is defined
    ::ForceFieldToolKit::SharedFcns::checkElementPDB

    # parse the geometry optimization output file
    set inFile [open $logFile r]

    while { ![eof $inFile] } {
        set inLine [string trim [gets $inFile]]
        if { $inLine eq "OPTKING Finished Execution" } {
            # burn the coord header
            for {set i 0} {$i < 8} {incr i} { gets $inFile }
            # read coordinates
            set coords {}
            while { ![regexp {^-*$} [set inLine [string trim [gets $inFile]]]] } {
                lappend coords [lrange $inLine 1 4]
            }
            set coords [::ForceFieldToolKit::SharedFcns::LonePair::addLPCoordinate $coords]

            # (re)set the coords
            set temp [atomselect $molId all]
            $temp set {x y z} $coords
            $temp delete
            unset coords
        } else {
            continue
        }
    }

    # write the new coords to file
    [atomselect $molId all] writepdb $optPdb

    # clean up
    close $inFile
    mol delete $molId

}
#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::auxFit_ChargeOpt { refmolid resName } {
     # Empty proc, used for other QM softwares (i.e. ORCA)
}

#===========================================================================================================
# CHARGE PARAMETRIZATION: WATER INTERACTION
#===========================================================================================================

proc ::ForceFieldToolKit::Psi4::genZmatrix {outFolderPath basename donList accList qmProc qmMem qmRoute qmCharge qmMult } {
    # NEW 02/01/2019: Checking "element" is defined
    ::ForceFieldToolKit::SharedFcns::checkElementPDB

    # computes the geometry-dependent position of water and
    # writes Psi4 input files for optimizing water interaction
    # assign Psi4 atom names and gather x,y,z for output python file
    set Gnames {}
    set atom_info {}
    for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
        set temp [atomselect top "index $i"]
        lappend atom_info [list [$temp get element] [$temp get x] [$temp get y] [$temp get z]]
        lappend Gnames [$temp get element][expr $i+1]
        $temp delete
    }

    # length of $atom_info
    set len_mol [llength $atom_info]

    #========#
    # DONORS #
    #========#
    foreach donorAtom $donList {
		set donorName [[atomselect top "index $donorAtom"] get name]

		# open output file
        set outname [file join $outFolderPath ${basename}-DON-${donorName}.py]
        set outfile [open $outname w]

        # start python code
        puts $outfile "import psi4"
        puts $outfile "import optking"
        puts $outfile "import qcelemental as qcel"
        puts $outfile "import numpy as np"
        puts $outfile ""
        puts $outfile {# Summary: optimizes intermolecular coordinates between two frozen monomers.}
        puts $outfile {# Split dimers in input with "--"}
        puts $outfile ""
        puts $outfile "psi4.set_num_threads($qmProc)"
        puts $outfile "psi4.set_output_file(\"${basename}-DON-${donorName}.out\", False)"
        puts $outfile {molecule = psi4.geometry("""}

		# write the cartesian coords for the molecule
		foreach atom_entry $atom_info {
		    puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
		}

        puts $outfile "--"

	    # build / write the zmatrix
        set ref_label [list A1 A2 A3 B1 B2 B3]
        set ref_val [::ForceFieldToolKit::Psi4::writeZmat $donorAtom donor $Gnames $outfile $len_mol]
        for {set i 0} {$i < 6} {incr i} {
            set [lindex $ref_label $i] [lindex $ref_val $i]
        }

        puts $outfile "nocom"
        puts $outfile "unit angstrom"
        puts $outfile {""")}
        puts $outfile ""

        ::ForceFieldToolKit::Psi4::write_optZmat $qmMem $qmMult $qmCharge $qmRoute $len_mol $outfile $A1 $A2 $A3 $B1 $B2 $B3

        close $outfile
    }; # end DONORS

    #===========#
    # ACCEPTORS #
    #===========#
    foreach acceptorAtom $accList {
		set acceptorName [[atomselect top "index $acceptorAtom"] get name]

		# open output file
		set outname [file join $outFolderPath ${basename}-ACC-${acceptorName}.py]
		set outfile [open $outname w]

        # start python code
        puts $outfile "import psi4"
        puts $outfile "import optking"
        puts $outfile "import qcelemental as qcel"
        puts $outfile "import numpy as np"
        puts $outfile ""
        puts $outfile {# Summary: optimizes intermolecular coordinates between two frozen monomers.}
        puts $outfile {# Split dimers in input with "--"}
        puts $outfile ""
        puts $outfile "psi4.set_num_threads($qmProc)"
        puts $outfile "psi4.set_output_file(\"${basename}-ACC-${acceptorName}.out\", False)"
        puts $outfile "molecule = psi4.geometry(\"\"\""

		# write the cartesian coords for the molecule
		foreach atom_entry $atom_info {
		    puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
		}

        puts $outfile "--"

		# build / write the zmatrix
        set ref_label [list A1 A2 A3 B1 B2 B3]
        set ref_val [::ForceFieldToolKit::Psi4::writeZmat $acceptorAtom acceptor $Gnames $outfile $len_mol]
        for {set i 0} {$i < 6} {incr i} {
            set [lindex $ref_label $i] [lindex $ref_val $i]
        }

        puts $outfile "nocom"
        puts $outfile "unit angstrom"
        puts $outfile {""")}
        puts $outfile ""

        ::ForceFieldToolKit::Psi4::write_optZmat $qmMem $qmMult $qmCharge $qmRoute $len_mol $outfile $A1 $A2 $A3 $B1 $B2 $B3

        close $outfile

		# CARBONYL AND HALOGEN EXCEPTIONS
		# there is a special exception for X=O cases (e.g. C=O, P=O, S=O)
		# where we need to write two additional files
        # another speical exception for halogens, writing one additional file
		::ForceFieldToolKit::GenZMatrix::writeExceptionZMats $acceptorName $acceptorAtom $Gnames $atom_info

    }
}
#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::writeZmat { aInd class gnames outfile len_mol } {
	# builds and writes the z-matrix to file

	# passed:
	#		aInd = atom index for interacting atom
	#		class = acceptor / donor
	#		gnames = list of Psi4-style names
	#		outfile = file handle where output is written

	# returns: nothing

	# make a selection for atom A (the interaction site), and find all attached atoms
	set aSel [atomselect top "index $aInd"]
	set bondlistA [lindex [$aSel getbonds] 0]
  # puts $aInd
  # puts $aSel
  # puts $bondlistA

	# z-matrix will be different for n=1 and n>1 (each will have their own special cases as well)
	if { [llength $bondlistA] == 1} {
		#=======#
		# N = 1 #
		#=======#
		# when n=1 there are two cases in which only 1D scanning is possible:
		# a diatomic molecule, or when C--B--A angle is 180 degrees (i.e. linear)

		set diatomic 0; # defaults to false
		set linear 1; # defaults to true

		# check for a diatomic molecule
		if { [molinfo top get numatoms] == 2} { set diatomic 1 }
		# check if C--B--A is linear
		if { !$diatomic } {
      set bInd [lindex $bondlistA 0]
			set bSel [atomselect top "index $bondlistA"]
			set bondlistB [lindex [$bSel getbonds] 0]
			foreach ele $bondlistB {
				# if C = A, skip
				if { $ele == $aInd } { continue }
				# check for a non-linear angle (+/- 2 degrees)
				if { [expr {abs([measure angle [list $aInd $bondlistA $ele]])}] <= 178.0 } {
					# found a non-linear C atom; unset the flag and stop looking
          set cInd $ele
					set linear 0; break
				} else {
					# keep looking
					continue
				}
			}
			# clean up atom selections
			$bSel delete
		}

    # define x and Ow as the indices (start from 1) of the dummy atom and O, respectively
    set x [expr $len_mol + 1]
    set Ow [expr $len_mol + 2]

		if { $diatomic || $linear } {
			# if either diatomic or linear, build a zmatrix for a 1D scan

			# get the relevant gnames
			set aGname [expr $aInd + 1]
			set bGname [expr $bondlistA + 1]
      set cGanme [expr $cInd + 1]   # added

			# for non-diatomic linear cases, we will define x in cartesian coords
			if { $linear } {
				set bSel [atomselect top "index $bondlistA"]
				set v1 [vecnorm [vecsub [measure center $aSel] [measure center $bSel]]]
				if { [lindex $v1 0] == 0 } {
					set v2 "1 0 0"
				} elseif { [lindex $v1 1] == 0 } {
					set v2 "0 1 0"
				} elseif { [lindex $v1 2] == 0 } {
					set v2 "0 0 1"
				} else {
					set v2 [list 1 [expr -1.0*[lindex $v1 0]/[lindex $v1 1]] 0]
				}
				set v2 [vecnorm $v2]
				set xPos [vecadd $v2 [measure center $aSel]]
			}

			# positioning of x is the same for donors and acceptors
			if { $diatomic } {
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s" x $aGname 1.0 $bGname 90.00]
			} else {
				puts $outfile [format "%3s  %.4f  %.4f  %.4f" x [lindex $xPos 0] [lindex $xPos 1] [lindex $xPos 2]]
			}

			# write the rest of the z-matrix
			if { $class eq "donor" } {
				# donor
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s %7s  %6s" O $aGname rAH $x 90.00 $bGname 180.00]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s %7s  %6s" H $Ow 0.9572 $aGname 127.74 $x 0.00]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s %7s  %6s\n" H $Ow 0.9572 $aGname 127.74 $x 180.00]
				puts $outfile "rAH = 2.0"
			} else {
				# acceptor
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s %7s  %6s" O $aGname rAH $x 90.00 $bGname 180.00]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s %7s  %6s" H $Ow 0.9572 $aGname 104.52 $x   0.00]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s %7s  %6s\n" H $Ow 0.9572 H2w 104.52 $x 0.0]
				puts $outfile "rAH = 2.0"
			}

			# clean up and return
			$aSel delete
			if { $linear } { $bSel delete }
      # assign reference atoms for diatomic/linear N = 1 case
      lassign "$aGname $bGname $cGname" A1 A2 A3
      lassign "[expr $len_mol + 1] [expr $len_mol + 2] [expr $len_mol + 3]" B1 B2 B3

			# done with special n=1 cases (diatomic or linear)
		} else {
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
			set aGname [expr $aInd + 1]
			set bGname [expr $bInd + 1]
			set cGname [expr $cInd + 1]
      # puts "aGname $aGname bGname $bGname cGname $cGname"

			if { $class eq "donor" } {
				# donor
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" x $aGname 1.0 $bGname 90.00 $cGname dih]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" O $aGname rAH $x 90.00 $bGname 180.00]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" H $Ow 0.9572 $aGname 127.74 $x 0.00]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s\n" H $Ow 0.9572 $aGname 127.74 $x 180.00]
				puts $outfile "rAH = 2.0"
				puts $outfile "dih = 0.0"

        # assign reference atoms for "normal" N = 1 case donor
        lassign "$aGname $bGname $cGname" A1 A2 A3
        lassign "[expr $len_mol + 1] [expr $len_mol + 2] [expr $len_mol + 3]" B1 B2 B3
			} else {
				# acceptor
				# call helper function to find probe position
				set probePos [::ForceFieldToolKit::Psi4::placeProbe $aSel]
				# note that Psi4 doesn't like 180 degree angles, so we have to be a little clever
				# make some measurements of probe position
				set mAng [::QMtool::bond_angle $probePos [measure center $aSel] [measure center $cSel]]
				set mDih [::QMtool::dihed_angle $probePos [measure center $aSel] [measure center $cSel] [measure center $bSel]]

        # need to redefine x and Ow because the dummy atom X is at index $len_mol + 2
        # and Ow is at index $len_mol + 3 instead
        # In addition, define H1w as the index of H1w
        set x [expr {$len_mol + 2}]
        set Ow [expr {$len_mol + 3}]
        set H1w [expr {$len_mol + 1}]

				puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" H1w $aGname rAH $cGname [format %3.2f $mAng] $bGname [format %3.2f $mDih]]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" x $H1w 1.0 $aGname 90.00 $cGname 0.00]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" Ow $H1w 0.9527 $x 90.00 $aGname 180.00]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s\n" H2w $Ow 0.9527 $H1w 104.52 $x dih]
				puts $outfile "rAH = 2.0"
				puts $outfile "dih = 0.0"

        # assign reference atoms for "normal" N = 1 case aceptor
        lassign "$aGname $bGname $cGname" A1 A2 A3
        lassign "[expr $len_mol + 1] [expr $len_mol + 2] [expr $len_mol + 3]" B1 B2 B3
			}

			# clean up and return
			$aSel delete; $bSel delete; $cSel delete
		}; # end of 'normal' N=1

	} else {
		#=======#
		# N > 1 #
		#=======#

		# find some information about atom B
		set bInd [lindex $bondlistA 0]
		set bSel [atomselect top "index $bInd"]
		set cInd [lindex $bondlistA 1]
		set cSel [atomselect top "index $cInd"]

		# find gnames
		set aGname [expr $aInd + 1]
		set bGname [expr $bInd + 1]
		# set cGname [expr $cInd + 1]   # if cGname is set to be "x2" below, then this line is redundant

		# test if C is valid choice
		set abcAng [expr {abs([measure angle [list $aInd $bInd $cInd]])}]
		if { $abcAng < 2 || $abcAng > 178 } { set validC 0 } else {	set validC 1 }
		unset abcAng

		# find probe coords
		set probePos [::ForceFieldToolKit::Psi4::placeProbe $aSel]
		set mAng [::QMtool::bond_angle $probePos [measure center $aSel] [measure center $bSel]]

		if { !$validC } {
			# if C is invalid, ABC are linear and we need a second dummy atom in lieu of the original C atom
			set cGname [expr $len_mol + 1]
			set x2Pos [coordtrans [trans center [measure center $aSel] axis [vecsub [measure center $cSel] [measure center $bSel]] 180.0 deg] $probePos]
			set mDih [::QMtool::dihed_angle $probePos [measure center $aSel] [measure center $bSel] $x2Pos]
			puts $outfile [format "%3s  %.4f  %.4f  %.4f" X [lindex $x2Pos 0] [lindex $x2Pos 1] [lindex $x2Pos 2]]

      # define x, Ow, H1w, and H2w
      set x [expr $len_mol + 3]
      set Ow [expr $len_mol + 2]
      set H1w [expr $len_mol + 4]
      set H2w [expr $len_mol + 5]

		} else {
			# C is valid, we can use it to define the dihedral
			set mDih [::QMtool::dihed_angle $probePos [measure center $aSel] [measure center $bSel] [measure center $cSel]]

      set cGname [expr $len_mol + 1]

      # define x, Ow, H1w, and H2w
      set x [expr $len_mol + 2]
      set Ow [expr $len_mol + 1]
      set H1w [expr $len_mol + 3]
      set H2w [expr $len_mol + 4]
		}

		if { $class eq "donor" } {
			# donor
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" Ow $aGname rAH $bGname [format %3.2f $mAng] $cGname [format %3.2f $mDih]]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" x $Ow 1.0 $aGname 90.00 $bGname dih]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" H1w $Ow 0.9572 $aGname 127.74 $x 0.00]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s\n" H2w $Ow 0.9572 $aGname 127.74 $x 180.00]
            puts $outfile "rAH = 2.0"
            puts $outfile "dih = 0.0"
		} else {
			# acceptor
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" H1w $aGname rAH $bGname [format %3.2f $mAng] $cGname [format %3.2f $mDih]]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" x $H1w 1.0 $aGname 90.00 $bGname 0.00]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" Ow $H1w 0.9572 $x 90.00 $aGname 180.00]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s\n" H2w $Ow 0.9572 $H1w 104.52 $x dih]
			puts $outfile "rAH = 2.0"
			puts $outfile "dih = 0.0"
		}

		# clean up atomselections and return
		$aSel delete; $bSel delete; $cSel delete

    # assign the reference atoms for N > 1 case
    lassign "$aGname $bGname $cGname" A1 A2 A3
    lassign "[expr $len_mol + 1] [expr $len_mol + 2] [expr $len_mol + 3]" B1 B2 B3
	}
  # return the reference atoms
  # A1, A2, A3 and B1, B2, B3 are defined in
  # https://github.com/zachglick/fftk-psi4-translation/blob/main/psi4/charge-opt-water/PRLD-DON-H1.py
  return [list $A1 $A2 $A3 $B1 $B2 $B3]
}
#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::write_optZmat {qmMem qmMult qmCharge qmRoute len_mol outfile A1 A2 A3 B1 B2 B3} {

  # follow https://github.com/zachglick/fftk-psi4-translation/blob/main/psi4/charge-opt-water
  puts $outfile {# Ordering of reference atoms is essential for properly defining contraints.}
  puts $outfile {# Note that if either of the theta angles are linear. optking will crash.}
  puts $outfile {# Keep to the ordering below and everything should be fine :)}
  puts $outfile {#  A1: Ligand donor/acceptor atom}
  puts $outfile {#  A2: Ligand atom in-plane with dih (NOT in-line with H-bond)}
  puts $outfile {#  A3: Ligand atom}
  puts $outfile {#  B1: Water donor/acceptor atom}
  puts $outfile {#  B2: Water atom in-plane with dih (NOT in-line with H-bond)}
  puts $outfile {#  B3: Final water atom}
  puts $outfile ""

  puts $outfile {# Here are the dimer coordinates that are used with their definitions.}
  puts $outfile {#  R         : Distance A1 to B1 (rAH)}
  puts $outfile {#  theta_A   : Angle,          A2-A1-B1}
  puts $outfile {#  theta_B   : Angle,          A1-B1-B2}
  puts $outfile {#  tau       : Dihedral angle, A2-A1-B1-B2 (dih)}
  puts $outfile {#  phi_A     : Dihedral angle, A3-A2-A1-B1}
  puts $outfile {#  phi_B     : Dihedral angle, A1-B1-B2-B3}
  puts $outfile ""

  puts $outfile {# Coords are 1-indexed without dummies here}
  puts $outfile "dimer = {"
  puts $outfile "    \"Natoms per frag\": \[$len_mol, 3\],"
  puts $outfile {    "A Frag": 1,}
  puts $outfile "    \"A Ref Atoms\": \[\[$A1\], \[$A2\], \[$A3\]\],"
  puts $outfile {    "A Label": "PRLD",}
  puts $outfile {    "B Frag": 2,}
  puts $outfile "    \"B Ref Atoms\": \[\[$B1\], \[$B2\], \[$B3\]\],"
  puts $outfile {    "B Label": "Water",}
  puts $outfile {    "Frozen": ["theta_A", "theta_B", "phi_A", "phi_B"],}
  puts $outfile "}"
  puts $outfile ""

  puts $outfile "psi4.set_memory(\"$qmMem GB\")"
  puts $outfile "molecule.set_multiplicity($qmMult)"
  puts $outfile "molecule.set_molecular_charge($qmCharge)"
  puts $outfile {psi4.core.clean_options()}
  puts $outfile "psi4_options = {"
  puts $outfile {    "basis": "6-31G*",}
  puts $outfile {    "d_convergence": 9,}
  puts $outfile {    "frag_mode": "multi",}
  puts $outfile {    "freeze_intrafrag": 'true',}
  puts $outfile "}"
  puts $outfile {psi4.set_options(psi4_options)}
  puts $outfile "result = optking.optimize_psi4(\"$qmRoute\", **{\"interfrag_coords\": str(dimer)})"
  puts $outfile ""

  puts $outfile {# Can use qcel to output final rAH and dih values}
  puts $outfile {xyzs = result["final_molecule"]["geometry"]}
  puts $outfile {xyzs = np.array(xyzs)}
  puts $outfile {xyzs = np.reshape(xyzs, (-1,3))}
  puts $outfile {# coords are zero-indexed here:}
  puts $outfile "rAH = qcel.util.measure_coordinates(xyzs, $A1, $B1, True) # in bohr"
  puts $outfile "# dih = qcel.util.measure_coordinates(xyzs, \[...\], True) # in degrees   # dih does not exist in every case"
  puts $outfile ""

  puts $outfile {bohr_to_ang = qcel.constants.conversion_factor("bohr", "angstrom")}
  puts $outfile {print(rAH * bohr_to_ang)}
  puts $outfile {print(dih)}
}
#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::placeProbe { aSel } {
	# helper function for genZmatrix
	# computes the position of the first probe atom
	# (i.e. Oxygen for HB donors or Hydrogen for HB acceptors)

	# takes in atomselection for the interacting atom
	# returns x,y,z of the first probe atom (i.e. O or H)

	# setup some variables
	set mol [$aSel molid]
	set aCoord [measure center $aSel]

	set bonded [lindex [$aSel getbonds] 0]
	set numBonded [llength $bonded]

	set all [atomselect $mol "all"]
	set allCenter [measure center $all]


	set count 0
	set normavg {0 0 0}


	# compute the probe position based on geometry
	# n = 1, i.e. b---a
	if { $numBonded == 1 } {
		set temp1 [atomselect $mol "index [lindex $bonded 0]"]
		set normavg [vecsub $aCoord [measure center $temp1]]
		incr count
		$temp1 delete
	}


	# n = 2, i.e. b1---a---b2
	if { $numBonded == 2 } {
		set temp1 [atomselect $mol "index [lindex $bonded 0]"]
		set bondvec1 [vecnorm [vecsub $aCoord [measure center $temp1]]]
		set temp2 [atomselect $mol "index [lindex $bonded 1]"]
		set bondvec2 [vecnorm [vecsub $aCoord [measure center $temp2]]]

		# check for linearity, project point away from rest of molecule
        if { abs([vecdot $bondvec1 $bondvec2]) > 0.95 } {
            # check that center of molecule doesn't already lie on our line, or worse, our atom
            # if it does, we don't care about where we place the water relative to the rest of
            # molecule
            set flag 0
            if { [veclength [vecsub $allCenter $aCoord]] < 0.01 } {
                set flag 1
            } elseif { abs([vecdot $bondvec1 [vecnorm [vecsub $allCenter $aCoord]]]) > 0.95 } {
                set flag 1
            }
            if { $flag } {
                if { [lindex $bondvec1 0] == 0 } {
                   set probeCoord] "1 0 0"
                } elseif { [lindex $bondvec1 1] == 0 } {
                   set probeCoord "0 1 0"
                } elseif { [lindex $bondvec1 2] == 0 } {
                   set probeCoord "0 0 1"
                } else {
                   set probeCoord [list 1 [expr -1.0*[lindex $bondvec1 0]/[lindex $bondvec1 1]] 0]
                }
                set probeCoord [vecadd $aCoord [vecscale 2.0 [vecnorm $probeCoord]]]
                return $probeCoord
            }
            set alpha [vecdot $bondvec1 [vecsub $allCenter $aCoord]]
            # same side as mol center
            set probeCoord [vecsub $allCenter [vecscale $alpha $bondvec1]]
            # reflect, extend
            set probeCoord [vecadd $aCoord [vecscale 2.0 [vecnorm [vecsub $aCoord $probeCoord]]]]
            return $probeCoord
        } else {
            set normavg [vecadd $bondvec1 $bondvec2]
        }

		incr count
		$temp1 delete; $temp2 delete; $all delete
	}

	# n > 2, there are many cases; majority are n=3 and n=4
	if { $numBonded > 2 } {
		# cycle through to find all normals
	    for {set i 0} {$i <= [expr $numBonded-3]} {incr i} {

	        set temp1 [atomselect $mol "index [lindex $bonded $i]"]
	        # normalize bond vectors first
	        set normPos1 [vecadd $aCoord [vecnorm [vecsub [measure center $temp1] $aCoord]]]

	        for {set j [expr $i+1]} {$j <= [expr $numBonded-2]} {incr j} {
				set temp2 [atomselect $mol "index [lindex $bonded $j]"]
				set normPos2 [vecadd $aCoord [vecnorm [vecsub [measure center $temp2] $aCoord]]]

				for {set k [expr $j+1]} {$k <= [expr $numBonded-1]} {incr k} {
					set temp3 [atomselect $mol "index [lindex $bonded $k]"]
					set normPos3 [vecadd $aCoord [vecnorm [vecsub [measure center $temp3] $aCoord]]]

					# get the normal vector to the plane formed by the three atoms
					set vec1 [vecnorm [vecsub $normPos1 $normPos2]]
					set vec2 [vecnorm [vecsub $normPos2 $normPos3]]
					set norm [veccross $vec1 $vec2]

					# check that the normal vector and atom of interest are on the same side of the plane
					set d [expr -1.0*[vecdot $norm $normPos1]]
					if { [expr $d + [vecdot $norm $aCoord]] < 0 } {
						set norm [veccross $vec2 $vec1]
					}

					# will average normal vectors at end
					set normavg [vecadd $normavg $norm]
					incr count
					$temp3 delete
				}
				$temp2 delete
	        }
	        $temp1 delete
	    }
	}

	# finish up and return
	set normavg [vecscale [expr 1.0/$count] $normavg]
	set probeCoord [vecadd $aCoord [vecscale 2.0 [vecnorm $normavg]]]
	return $probeCoord
}
#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::write120filesWI { outFolderPath basename aName aInd atom_info aGname bGname cGname qmProc qmMem qmRoute qmCharge qmMult } {
    # writes single point energy files required for charge optimization
    # hard coded for HF/6-31G* and MP2/6-31G*

    # write two slightly different files
    foreach altPos {"a" "b"} dihed {180 0} {

    	# open output file
    	set outname [file join $outFolderPath ${basename}-ACC-${aName}-120${altPos}.gau]
    	set outfile [open $outname w]

    	# write the header
    	puts $outfile "%chk=${basename}-ACC-${aName}-120${altPos}.chk"
    	puts $outfile "%nproc=$qmProc"
    	puts $outfile "%mem=${qmMem}GB"
      puts $outfile "$qmRoute"
    	puts $outfile ""
    	puts $outfile "<qmtool> simtype=\"Geometry optimization\" </qmtool>"
    	puts $outfile "${basename}-ACC-${aName}-120${altPos}"
    	puts $outfile ""
    	puts $outfile "$qmCharge $qmMult"

    	# write the cartesian coords for the molecule
    	foreach atom_entry $atom_info {
    	    puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
    	}

    	# write custom zmatrix
    	set ang 120

    	puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" H1w $aGname rAH $bGname [format %3.2f $ang] $cGname [format %3.2f $dihed]]
    	puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" x H1w 1.0 $aGname 90.00 $bGname 0.00]
    	puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" Ow H1w 0.9527 x 90.00 $aGname 180.00]
    	puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s\n" H2w Ow 0.9527 H1w 104.52 x dih]
    	puts $outfile "rAH  2.0"
    	puts $outfile "dih  0.0"

    	# close up
    	close $outfile
    }
}
#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::write90filesWI { outFolderPath basename aName aInd atom_info aGname bGname cGname qmProc qmMem qmRoute qmCharge qmMult } {
    # writes single point energy files required for charge optimization
    # hard coded for HF/6-31G* and MP2/6-31G*

    # open output file
    set outname [file join $outFolderPath ${basename}-ACC-${aName}-ppn.gau]
    set outfile [open $outname w]

    # write the header
    puts $outfile "%chk=${basename}-ACC-${aName}-ppn.chk"
    puts $outfile "%nproc=$qmProc"
    puts $outfile "%mem=${qmMem}GB"
  puts $outfile "$qmRoute"
    puts $outfile ""
    puts $outfile "<qmtool> simtype=\"Geometry optimization\" </qmtool>"
    puts $outfile "${basename}-ACC-${aName}-ppn"
    puts $outfile ""
    puts $outfile "$qmCharge $qmMult"

    # write the cartesian coords for the molecule
    foreach atom_entry $atom_info {
        puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
    }

    # write custom zmatrix
    puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" H1w $aGname rAH $bGname 90.0 $cGname 90.0]
    puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" x H1w 1.0 $aGname 90.00 $bGname 0.00]
    puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" Ow H1w 0.9527 x 90.00 $aGname 180.00]
    puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s\n" H2w Ow 0.9527 H1w 104.52 x dih]
    puts $outfile "rAH  2.0"
    puts $outfile "dih  0.0"

    # close up
    close $outfile
}
#===========================================================================================================

proc ::ForceFieldToolKit::Psi4::writeSPfilesWI { outFolderPath basename qmProc qmMem qmCharge qmMult } {
    # writes single point energy files required for charge optimization
    # hard coded for HF/6-31G* and MP2/6-31G*

    for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
        set temp [atomselect top "index $i"]
        lappend atom_info [list [$temp get element] [$temp get x] [$temp get y] [$temp get z]]
        lappend Gnames [$temp get element][expr $i+1]
        $temp delete
    }
    mol delete top

    ### Write the HF Single Point File
    # open output file
    set outfile [open [file join $outFolderPath ${basename}-sp-HF.py] w]

    # start python code
    puts $outfile "import psi4"
    puts $outfile ""
    puts $outfile "psi4.set_output_file(\"${basename}-sp-HF.out\", False)"
    puts $outfile "psi4.set_memory(\"$qmMem GB\")"
    puts $outfile ""
    puts $outfile {molecule = psi4.geometry("""}

    # write the coordinates
    foreach atom_entry $atom_info {
        puts $outfile "[lindex $atom_entry 0] [format %16.8f [lindex $atom_entry 1]] [format %16.8f [lindex $atom_entry 2]] [format %16.8f [lindex $atom_entry 3]]"
   }

    puts $outfile {""")}
    puts $outfile ""
    puts $outfile "molecule.set_multiplicity($qmMult)"
    puts $outfile "molecule.set_molecular_charge($qmCharge)"
    puts $outfile "psi4.optimize(\"HF/6-31G*\", molecule=molecule)"

    puts $outfile ""
    close $outfile

    ### Write the MP2 Single Point File
    # open output file
    set outfile [open [file join $outFolderPath ${basename}-sp-MP2.py] w]

    puts $outfile "import psi4"
    puts $outfile ""
    puts $outfile "psi4.set_output_file(\"${basename}-sp-MP2.out\", False)"
    puts $outfile "psi4.set_memory(\"$qmMem GB\")"
    puts $outfile ""
    puts $outfile {molecule = psi4.geometry("""}

    # write the coordinates
    foreach atom_entry $atom_info {
    puts $outfile "[lindex $atom_entry 0] [format %16.8f [lindex $atom_entry 1]] [format %16.8f [lindex $atom_entry 2]] [format %16.8f [lindex $atom_entry 3]]"
    }

    puts $outfile {""")}
    puts $outfile ""
    puts $outfile "molecule.set_multiplicity($qmMult)"
    puts $outfile "molecule.set_molecular_charge($qmCharge)"
    puts $outfile "psi4.optimize(\"MP2/6-31G*\", molecule=molecule)"

    puts $outfile ""
    close $outfile

    ### Write TIP3P Single Point File
    set outfile [open [file join $outFolderPath wat-sp.py] w]
    puts $outfile "import psi4"
    puts $outfile ""
    puts $outfile "psi4.set_memory(\"$qmMem GB\")"
    puts $outfile "# RHF/6-31G* SCF=Tight"
    puts $outfile ""
    puts $outfile "# wat-sp"
    puts $outfile ""
    puts $outfile {molecule = psi4.geometry("""}
    puts $outfile "O"
    puts $outfile "H 1 0.9572"
    puts $outfile "H 1 0.9572 2 104.52"
    puts $outfile {""")}
    puts $outfile ""
    puts $outfile "molecule.set_multiplicity($qmMult)"
    puts $outfile "molecule.set_molecular_charge($qmCharge)"
    puts $outfile "psi4.optimize(\"RHF/6-31G*\", molecule=molecule)"
    close $outfile

}
#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::loadCOMFile { comfile } {
  # New QM Input loader
  set inFile [open $comfile r]
  set tmpFile [open "$comfile-tmp" w]

  set read_coords 0
  while { ![eof $inFile] } {
      set inLine [string trim [gets $inFile]]

      if { [string match {*psi4.geometry*} $inLine] } {
          set read_coords 1
          continue
      }

      if {$read_coords} {
          if { [string match {""")} $inLine] } {
              break
          }
          if {[llength $inLine] >= 4
              && [string is alpha [string index [lindex $inLine 0] 0] ]
              && [string is double [lindex $inLine 1]]} {
              set inLine [regsub "rAH" $inLine "2.0"]
              set inLine [regsub "dih" $inLine "0.0"]
              puts $tmpFile $inLine
          }
      }
  }
  close $inFile
  close $tmpFile
  ::QMtool::read_zmtfile "$comfile-tmp"
  file delete "$comfile-tmp"
  set molId [molinfo top]
  return $molId
}
#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::loadLOGFile { logfile } {
    # New QM Output loader
    set inFile [open $logfile r]
    # Create a temporary xyz file
    set firstXYZ 1
    while { ![eof $inFile] } {
        set inLine [string trim [gets $inFile]]
        if { [string match {==> Geometry <==} $inLine] } {
            # jump to the coordinates
            for {set i 0} {$i < 8} {incr i} {  # burn-in the header
                gets $inFile
            }
            # read coordinates
            set coords {}
            # get ligand coord
            while { [regexp {[A-Z]} [set inLine [string trim [gets $inFile]]]] } {
                lappend coords [lrange $inLine 0 3]
            }
            if { $firstXYZ eq 1 } { ;# for the first geom opt iteration, create xyz, fill and load it in VMD
                set tempFile [open temp.xyz w] ;# create the temp.xyz
                set numAtoms [llength $coords]
                puts $tempFile $numAtoms
                puts $tempFile ""
                for {set i 0} {$i < $numAtoms} {incr i} {
                     puts $tempFile [lindex $coords $i]
                }
                close $tempFile
                set molId [mol new temp.xyz] ;# load the temporary xyz file in VMD
                set firstXYZ 0
            } elseif { $firstXYZ eq 0 } { ;# for the following opt iterations addfiles
                # add a new frame, set the coords
                mol addfile temp.xyz
                for {set i 0} {$i < [llength $coords]} {incr i} {
                    set temp [atomselect $molId "index $i"]
                    $temp set x [lindex $coords $i 1]
                    $temp set y [lindex $coords $i 2]
                    $temp set z [lindex $coords $i 3]
                    $temp delete
                }
            }
            unset coords
        }
    }
    file delete temp.xyz ;# remove the temporary xyz file

    close $inFile
    return $molId
}

#===========================================================================================================
#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::getscf_ChargeOpt { file simtype } {
   # Get SCF energies from Psi4 Output
   set scfenergies {}

   set fid [open $file r]

   set hart_kcal 1.041308e-21; # hartree in kcal
   set mol 6.02214e23;

   set num 0
   set ori 0
   set tmpscf {}
   set optstep 0
   set scanpoint 0

   while {![eof $fid]} {
     set line [string trim [gets $fid]]

     if {[string match "Total Energy =" $line]} {
         set scf [lindex $line 2]
         set scfkcal [expr {$scf*$hart_kcal*$mol}]
         if {$num==0} { set ori $scf }
         set scfkcalori [expr {($scf-$ori)*$hart_kcal*$mol}]
         set tmpscf [list $num $scfkcal]

         incr num
     }

   }
   close $fid
   if {[llength $tmpscf]} { lappend scfenergies $tmpscf }

   return $scfenergies
}
#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::getMolCoords_ChargeOpt { file numMolAtoms } {

   set fid [open $file r]
   while { ![eof $fid] } {
      set inLine [string trim [gets $fid]]
      set coordlist {}   # define inside the while loop so they can be overwirtten until the last optimization
      if { $inLine eq "==> Geometry <==" } {
         # jump to coordinates
         for {set i 0} {$i < 8} {incr i} {  # burn-in the header
             gets $inFile
         }
         # read coordinates
         while { [regexp {[A-Z]} [set inLine [string trim [gets $fid]]]] } {
            lappend coordlist [lrange $inLine 1 3]
         }
      } else {
          continue
      }
  }

  close $fid

  return [lrange $coordlist 0 [expr $numMolAtoms - 1]]

}
#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::getWatCoords_ChargeOpt { file } {

   set fid [open $file r]
   while { ![eof $fid] } {
      set inLine [string trim [gets $fid]]
      set coordlist {}   # define inside the while loop so they can be overwirtten until the last optimization
      set atomnames {}
      if { $inLine eq "==> Geometry <==" } {
         # jump to coordinates
         for {set i 0} {$i < 8} {incr i} {  # burn-in the header
             gets $inFile
         }
         # read coordinates
         while { [regexp {[A-Z]} [set inLine [string trim [gets $fid]]]] } {
            lappend atomnames [lindex $inLine 0]
            lappend coordlist [lrange $inLine 1 3]
         }
      } else {
          continue
      }
  }
  set numAtoms [llength $atomnames]
  set Hcount 0
  for {set i [expr $numAtoms - 3]} {$i < $numAtoms} {incr i} {
      set name [lindex $atomnames $i]
      if { [string match "O*" $name] } {
         set Ocoord [lindex $coordlist $i]
      } elseif { [string match "H*" $name] && $Hcount == 1} {
         set H2coord [lindex $coordlist $i]
         set Hcount 2
      } elseif { [string match "H*" $name] && $Hcount == 0} {
         set H1coord [lindex $coordlist $i]
         set Hcount 1
      }
  }

  close $fid

  set coords [list $Ocoord $H1coord $H2coord]
  return $coords
}
#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::getDipoleData_ChargeOpt { filename } {

   set fid [open $file r]
   while { ![eof $fid] } {
      set inLine [string trim [gets $fid]]

      # define inside the while loop so they can be overwirtten until the last optimization
      set coor {}
      set qmVec {}
      set qmMag {}

      if { $inLine eq "==> Geometry <==" } {
         # jump to coordinates
         for {set i 0} {$i < 8} {incr i} {  # burn-in the header
             gets $inFile
         }
         # read coordinates
         while { [regexp {[A-Z]} [set inLine [string trim [gets $fid]]]] } {
            lappend coor [lrange $inLine 1 3]
         }
       } elseif { $inLine eq "Dipole Moment: [D]" } {
           set inLine [string trim [gets $inFile]]
           set qmVec [list [lindex $inLine 1] [lindex $inLine 3] [lindex $inLine 5]]
           set qmMag [lindex $inLine 7]
       } else {
           continue
       }
   }

   # we're done with the output file
   unset inLine
   close $fid

   return [list $coor $qmVec $qmMag]
}

#======================================================
proc ::ForceFieldToolKit::Psi4::resetDefaultsGenZMatrix {} {
    # resets Psi4 settings to default for water interaction tab

    # reset
    set ::ForceFieldToolKit::GenZMatrix::qmProc 1
    set ::ForceFieldToolKit::GenZMatrix::qmMem 1
    set ::ForceFieldToolKit::GenZMatrix::qmCharge 0
    set ::ForceFieldToolKit::GenZMatrix::qmMult 1
    set ::ForceFieldToolKit::GenZMatrix::qmRoute "mp2/6-31G*"

}
#======================================================


#===========================================================================================================
# CHARGE PARAMETRIZATION: ELECTROSTATIC POTENTIAL
#===========================================================================================================

proc ::ForceFieldToolKit::Psi4::writeGauFile_ESP { chk qmProc qmMem qmCharge qmMult qmRoute gau } {

    ### make sure that the chk file is from Psi4 or that a pdb file was provided ###
    # Give error message if no suitable file was found
    if { [file extension $chk] ne ".chk" && [file extension $chk] ne ".pdb" } {
        tk_messageBox -type ok -icon warning -message "Action halted on error!" -detail "The structure file was not recognized.\
                 Please provide a Psi4 Checkpoint file (.chk) or a PDB file (.pdb)."
        return
    }

    ### Get atomic coordinates ###
    if { [file extension $chk] eq ".pdb" } {
    # if a pdb file was provided extract the atomic information from there.
        mol new $chk
        # NEW 02/01/2019: Checking "element" is defined
        ::ForceFieldToolKit::SharedFcns::checkElementPDB
        # Get coordinates from VMD
        set atom_info {}
        for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
            set temp [atomselect top "index $i"]
            lappend atom_info [list [$temp get element] [$temp get x] [$temp get y] [$temp get z]]
            $temp delete
        }
        # remove loaded molecule
        mol delete top
        # write the .gau file
	set newCHKname "[file rootname $gau].chk"
        set gauFile [open $gau w]
        # write header
        puts $gauFile "%chk=[file tail $newCHKname]"
        puts $gauFile "%nproc=$qmProc"
        puts $gauFile "%mem=${qmMem}GB"
	# Modify route section to remove Geom=Checkpoint, which is not needed if pdb is used. Give message in tk.
	regsub -nocase {Geom=Checkpoint} $qmRoute {} qmRoute_mod
   	puts ""
        puts $gauFile "$qmRoute_mod"
        puts $gauFile ""
        puts $gauFile "<qmtool> simtype=\"ESP Calculation\" </qmtool>"
        puts $gauFile ""
        puts $gauFile "$qmCharge $qmMult"
        # write the coordinates
        foreach atom_entry $atom_info {
            puts $gauFile "[lindex $atom_entry 0] [format %16.8f [lindex $atom_entry 1]] [format %16.8f [lindex $atom_entry 2]] [format %16.8f [lindex $atom_entry 3]]"
        }
        # empty line to terminate
        puts $gauFile ""
        close $gauFile

    } else {
    # otherwise copy structural information from Psi4 checkpoint file
        # make a copy of the CHK file to prevent from overwriting the original
        set newCHKname "[file rootname $gau].chk"
        file copy -force $chk $newCHKname

         set gauFile [open $gau w]

         # write the .gau file
         puts $gauFile "%chk=[file tail $newCHKname]"
         puts $gauFile "%nproc=$qmProc"
         puts $gauFile "%mem=${qmMem}GB"
         puts $gauFile "$qmRoute"
         puts $gauFile ""
         puts $gauFile "<qmtool> simtype=\"ESP Calculation\" </qmtool>"
         puts $gauFile ""
         puts $gauFile "$qmCharge $qmMult"

         close $gauFile

    }

}
#===========================================================================================================

proc ::ForceFieldToolKit::Psi4::writeDatFile_ESP { inputName gauLog } {
        set count 1
        set auScale 0.52917720

        # make sure qmSoft variable is set to the right value for the selected output file
        if {[::ForceFieldToolKit::SharedFcns::checkWhichQM $gauLog]} {return}

        # name the file based on the inputName
        set datName ${inputName}.dat
        set datFile [open $datName w]

        set logFile [open $gauLog r]

	# read the number of atoms
	while { [lindex [set line [string map { \" {} } [gets $logFile]] ] 0] ne "NAtoms=" } {
		continue
        }

	set nAtoms [lindex $line 1]

	# read the number of fit centers
	while { [lrange [set line  [gets $logFile]] 1 8] ne "points will be used for fitting atomic charges" } {
		continue
        }

	set nFitCenters [lindex $line 0]

	puts $datFile "  $nAtoms  $nFitCenters"

	# go back to the beginning of the .log file
	seek $logFile 0

	# read the Atom Fit Centers
	while { [lrange [set line [string map { \" {} } [gets $logFile]]] 0 1] ne "Atomic Center" } {
		continue
        }
   set formatStr "                   %s   %s   %s"
#	set temp [format $formatStr [lindex $line 5] [lindex $line 6] [lindex $line 7]]
   puts $datFile [format $formatStr [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 5] / $auScale]] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 6] / $auScale]] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 7] / $auScale]]]
   while { [lrange [set line [gets $logFile]] 0 1] eq "Atomic Center" } {
       puts $datFile [format $formatStr [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 5] / $auScale]] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 6] / $auScale]] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 7] / $auScale]]]
		continue
    }

	# go back to the beginning of the .log file
	seek $logFile 0

	# read the fit values
	while { [lrange [set line [string map { \" {} } [gets $logFile]]] 0 3] ne "Electrostatic Properties (Atomic Units)" } {
		continue
    }
	# G09 Revision B does not give the fit values, so there must be a way to check this at this point of the code and give an error.
	# Otherwise the code will loop forever.
	while { [lindex [set line [gets $logFile]] 1] ne "Fit" && ![eof $logFile] } {
		continue
    }
        if { [eof $logFile] } {
            tk_messageBox -type ok -icon warning -message "Action halted on error!" -detail "Cannot find fitting values in log file. This may happen because you are using g09 revision B. Try to use a later version."
            return
        }
	lappend fitValue [lindex $line 2]
	while { [lindex [set line [gets $logFile]] 1] eq "Fit" } {
		lappend fitValue [lindex $line 2]
		continue
    }

	# go back to the beginning of the .log file
	seek $logFile 0

	# write formatted fit values to the data file
	while { [lrange [set line [string map { \" {} } [gets $logFile]]] 0 2] ne "ESP Fit Center" } {
		continue
   }
	set formatStr "   %s   %s   %s"
	puts -nonewline $datFile [format "   %s" [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [lindex $fitValue 0]]]
	puts $datFile [format $formatStr [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 6] / $auScale]] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 7] / $auScale]] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 8] / $auScale]]]
	while { [lrange [set line [gets $logFile]] 0 2] eq "ESP Fit Center" } {
		puts -nonewline $datFile [format "   %s" [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [lindex $fitValue $count]]]
		puts $datFile [format $formatStr [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 6] / $auScale]] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 7] / $auScale]] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 8] / $auScale]]]
		incr count
		continue
	}
	close $logFile
	close $datFile
}
#======================================================
proc ::ForceFieldToolKit::Psi4::resetDefaultsESP {} {

    # resets the Psi4 Settings to default values
    set ::ForceFieldToolKit::ChargeOpt::ESP::qmProc 1
    set ::ForceFieldToolKit::ChargeOpt::ESP::qmMem 1
    set ::ForceFieldToolKit::ChargeOpt::ESP::qmCharge 0
    set ::ForceFieldToolKit::ChargeOpt::ESP::qmMult 1
    set ::ForceFieldToolKit::ChargeOpt::ESP::qmRoute "#P HF/6-31G* SCF=Tight Geom=Checkpoint Pop=MK IOp(6/33=2,6/41=10,6/42=17)"

}

#===========================================================================================================
# BOND/ANGLE PARAMETRIZATION
#===========================================================================================================

proc ::ForceFieldToolKit::Psi4::WriteComFile_GenBonded { geomCHK com qmProc qmMem qmRoute qmCharge qmMult psf pdb lbThresh } {
    # Check if a file with .chk extension was given
    set fext [file extension $geomCHK]
    if { $fext ne ".chk" } {
        tk_messageBox -type ok -icon warning -message "Action halted on error!" -detail "File has not a \".chk\" extension. It might not be a Psi4 Checkpoint File."
        return
    }

    # load the molecule psf/pdb to get the internal coordinates
    set logID [::ForceFieldToolKit::SharedFcns::LonePair::loadMolExcludeLP $psf $pdb]
    ::QMtool::use_vmd_molecule $logID
    set zmat [::QMtool::modredundant_zmat]

    # make a copy of the CHK file to prevent Psi4 from overwriting the original
    set newCHKname "[file rootname $com].chk"
    file copy $geomCHK $newCHKname

    # assign atom names and gather x,y,z for the output file
    set atom_info {}
    for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
        set temp [atomselect top "index $i"]
        lappend atom_info [list [$temp get element] [$temp get x] [$temp get y] [$temp get z]]
        $temp delete
    }

    # write the com file
    set outfile [open $com w]

    # shamelessly stolen from qmtool
    set bonds {}
    set angles {}
    set dihedrals {}
    set lbList {}
    foreach entry $zmat {
        set indexes [lindex $entry 2]
        set type [string toupper [string index [lindex $entry 1] 0]]

        # check for linear angle
        if { $type eq "A" && [lindex $entry 3] > $lbThresh } {
            # angle qualifies as a "linear bend"
            set type "L"
            lappend lbList $indexes
        }

        # check if standard dihedrals are part of linear bend (undefined)
        set skipflag 0
        if { $type eq "D" && [llength $lbList] > 0 } {
            # test each linear bend in lbList against current dih indices definition
            foreach ang $lbList {
                # test forward and rev angle definitions
                if { [string match "*$ang*" $indexes] || [string match "*[lreverse $ang]*" $indexes] } {
                    # positive test -> leave this dihedral out
                    set skipflag 1
                    break
                }
            }
        }
        if { $skipflag } {
            continue
        }

        if {$type == "B"} {
            lappend bonds "([lindex $indexes 0], [lindex $indexes 1]),"
        } elseif { $type=="A" } {
            lappend angles "([lindex $indexes 0], [lindex $indexes 1], [lindex $indexes 2]),"
        } elseif { $type=="D" || $type=="I" || $type=="O" } {
            # impropers modeled as dihedrals because Gaussian ignores out-of-plane bends
            # assume Psi4 does the same thing
            lappend dihedrals "([lindex $indexes 0], [lindex $indexes 1], [lindex $indexes 2], [lindex $indexes 3]),"
        }

        # puts $outfile "$type $indexes A [regsub {[QCRM]} [lindex $entry 5] {}]"
    }

    puts $outfile "import numpy as np"
    puts $outfile "import psi4"
    puts $outfile "import optking"
    puts $outfile ""
    puts $outfile "psi4.core.set_output_file(\'hess.out\')"
    puts $outfile "psi4.set_memory(\'${qmMem}GB\')"
    puts $outfile "psi4.set_num_threads($qmProc)"
    puts $outfile ""

    puts $outfile "# Specify geometry here"
    puts $outfile "mol = psi4.geometry(\"\"\""
    puts $outfile "$qmCharge $qmMult"
    foreach atom_entry $atom_info {
        puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
    }
    puts $outfile "\"\"\")"
    puts $outfile ""

    # write the entry to the input file
    puts $outfile "# Specify internal coordinates  here"
    puts $outfile "# Note: atoms indices start at zero"
    puts $outfile "bonds = \["
    foreach bond $bonds {
        puts $outfile $bond
    }
    puts $outfile "]"
    puts $outfile ""

    puts $outfile "angles = \["
    foreach angle $angles {
        puts $outfile $angle
    }
    puts $outfile "]"
    puts $outfile ""

    puts $outfile "dihedrals = \["
    foreach dihedral $dihedrals {
        puts $outfile $dihedral
    }
    puts $outfile "]"
    puts $outfile ""

    puts $outfile "psi4_options  = {"
    puts $outfile "\"mp2_type\" : \"df\","
    puts $outfile "\"freeze_core\" : True,"
    puts $outfile "\"basis\" : \"6-31G*\","
    puts $outfile "\"g_convergence\" : \"gau_tight\","
    puts $outfile "}"
    puts $outfile ""

    puts $outfile "psi4.set_options(psi4_options)"
    puts $outfile ""

    puts $outfile "# Optimize and extract information, inc. final geometry"
    puts $outfile "json_output = optking.optimize_psi4(\'MP2\')"
    puts $outfile "E = json_output\[\'energies\']\[-1]"
    puts $outfile "print(f\"Optimized Energy: {E}\")"
    puts $outfile "xyz_array = np.array(json_output\[\'final_molecule\']\[\'geometry\'])"
    puts $outfile "xyz = xyz_array.reshape(mol.natom(), 3)"
    puts $outfile "mol.set_geometry(psi4.core.Matrix.from_array(xyz))"
    puts $outfile ""

    puts $outfile "optking.optparams.Params = optking.optparams.OptParams({})"
    puts $outfile ""

    puts $outfile "from optking import stre, bend, tors"
    puts $outfile "bonds = \[stre.Stre(*bond) for bond in bonds]"
    puts $outfile "angles = \[bend.Bend(*angle) for angle in angles]"
    puts $outfile "dihedrals = \[tors.Tors(*dihedral) for dihedral in dihedrals]"
    puts $outfile "coords = bonds + angles + dihedrals"
    puts $outfile ""

    puts $outfile "psi4.set_options(psi4_options)"
    puts $outfile "Z = \[mol.Z(i) for i in range(0,mol.natom())]"
    puts $outfile "masses = \[mol.mass(i) for i in range(0,mol.natom())]"
    puts $outfile "f1 = optking.frag.Frag(Z, xyz, masses, intcos=coords, frozen=False)"
    puts $outfile "OptMol = optking.molsys.Molsys(\[f1])"
    puts $outfile "print(OptMol)"
    puts $outfile ""

    puts $outfile "# Compute the Cartesian Hessian with psi4"
    puts $outfile "Hxyz = psi4.hessian(\'MP2\') # returns a psi4 matrix"
    puts $outfile "print(\'Cartesian hessian\')"
    puts $outfile "# print(Hxyz.to_array())"
    puts $outfile ""

    puts $outfile "# Transform hessian to internals with optking, returns a numpy array"
    puts $outfile "# Hq = optking.intcosMisc.hessian_to_internals(Hxyz.to_array(), OptMol)"
    puts $outfile "# print(\'Internal Coordinate Hessian\')"
    puts $outfile "# print(Hq)"
    puts $outfile ""

    puts $outfile "# Also print transformed Hessian to output File"
    puts $outfile "# psi4.core.Matrix.from_array(Hq).print_out()"
    puts $outfile ""

    puts $outfile "# Open the hess.out and append the bond/angle/dihedral indices and values to the end."
    puts $outfile "# They are needed for the Opt. Bonded tab"
    puts $outfile "with open('hess.out', 'a') as file:"
    puts $outfile "    file.write(\"***Results for Opt. Bonded***\\n\\n\")"
    puts $outfile "    file.write(\"{:5} {:15} {:20} \\n\".format('type', 'indices', 'ANG(DEG)'))"
    puts $outfile "    for frag in OptMol.fragments:"
    puts $outfile "        for i in range(len(frag.intcos)):"
    puts $outfile "            indices = str(frag.intcos\[i])\[3:-1].split(',')   # the indices goes like ' R(1,2)', ' D(1,2,3,4)', etc." 
    puts $outfile "            indices_zero = ''"
    puts $outfile "            for idx in indices:"
    puts $outfile "                indices_zero += str(int(idx) - 1) + ','"
    puts $outfile "            file.write(\"{:5} {:15} {:20} \\n\".format(str(frag.intcos\[i])\[1], '(' + indices_zero\[:-1] + ')', str(frag.q_show()\[i])))"
    puts $outfile "    file.write(\"\\n\")"
    puts $outfile "    file.write(\"***get_inthessian_kcal_BondAngleOpt***\\n\\n\")"
    puts $outfile "    Hq = OptMol.hessian_to_internals(Hxyz.to_array())"
    puts $outfile "    np.savetxt(file, Hq)"

    close $outfile

}

#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::zmatqm_BondAngleOpt { debug debugLog hessLogID hessLog {asktypeList 0} } {

    # localize forcefieldtoolkit debugging variables
    # variable debug

    # load the Psi4 Log files from the hessian calculation
    if { $debug } { puts -nonewline $debugLog "loading hessian log file..."; flush $debugLog }
    # set hessLogID [mol new $psf]

    # set molId [::ForceFieldToolKit::Psi4::loadLOGFile $logfile]

    set natoms [molinfo $hessLogID get numatoms]
    set nbonds 0
    set nangles 0
    set ndiheds 0
    set nimprops 0
    set zmatqm {0}
    set flag "Q"
    set scan {}
    set readOptParam 0
    set typeList {} ;# only for the get_inthessian_kcal proc

    set inFile [open "[file rootname $hessLog].out"]
    while { ![eof $inFile] } {
        set inLine [string trim [gets $inFile]]
        if { $inLine eq "***Results for Opt. Bonded***" } {
            # jump to the coordinates
            for { set i 0 } { $i < 2 } { incr i } { gets $inFile }
            while { [regexp {[0-9]} [set inLine [string trim [gets $inFile]]]] } {
                 set slist [concat {*}[split $inLine " ,()"]]
                 set typePar [lindex $slist 0]
                 lappend typeList $typePar
                 switch $typePar {
                   R { set typePar bond; incr nbonds }
                   B { set typePar angle; incr nangles }
                   D { set typePar dihed; incr ndiheds }
                   # What labels other than R, B, D does Psi4 have?
                   # T { set typePar dihed; incr ndiheds }
                   # O { set typePar imprp; incr nimprops }
                   # X { continue }
                   # Y { continue }
                   # Z { continue }
                 }
                 if { [string match "*bond" $typePar] } {
                     set indexlist [list [lindex $slist 1] [lindex $slist 2]]
                     set val [lindex $slist end]
                     lappend zmatqm [list "R$nbonds" $typePar $indexlist $val {{} {}} $flag $scan]

                 } elseif { [string equal "angle" $typePar] } {
                     set indexlist [list [lindex $slist 1] [lindex $slist 2] [lindex $slist 3]]
                     set val [lindex $slist end]
                     lappend zmatqm [list "A$nangles" $typePar $indexlist $val {{} {} {} {}} $flag $scan]

                 } elseif {[string equal "dihed" $typePar]} {
                     set indexlist [list [lindex $slist 1] [lindex $slist 2] [lindex $slist 3] [lindex $slist 4]]
                     set val [lindex $slist end]
                     lappend zmatqm [list "D$ndiheds" $typePar $indexlist $val {{} {} {}} $flag $scan]
                 }
            }
            set ncoords [expr { $nbonds + $nangles + $ndiheds + $nimprops }]
            set havepar 1
            set havefc 1
            lset zmatqm 0 [list $natoms $ncoords $nbonds $nangles $ndiheds $nimprops $havepar $havefc]
            #foreach ele $zmatqm { puts $ele }
        }
    }

    if { $asktypeList } {
        return $typeList ;# info required for get_inthessian_kcal proc
    } else {
        return $zmatqm ;# normal return for the proc
        #return [list $zmatqm $hessLogID]
    }

}
#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::get_inthessian_kcal_BondAngleOpt { {hessLogID ""} {hessLog ""} } {

    set fid [open "[file rootname $hessLog].out" r]
    set dimHess 0
    set hess {}
    set inthessian_kcal {}

    while {![eof $fid]} {
        set line [gets $fid]
        # jump to $hessian line
        if {[string first "***get_inthessian_kcal_BondAngleOpt***" $line] >= 0} {
            gets $fid
            while {![eof $fid]} {
                set rowdata [gets $fid]
                lappend hess $rowdata
            }
        }

    }

    # calling the zmatqm proc to obtain the typeList for the appropiate dim values
    set debug 0
    set debugLog ""
    set asktypeList 1
    set typeList [ zmatqm_BondAngleOpt $debug $debugLog $hessLogID $hessLog $asktypeList ]

    # convert hess from hartree*bohr2 to kcal*A2
    for { set i 0 } { $i < [llength $hess] } { incr i } {
        set rowdata [lindex $hess $i]
        set rowdata_kcal {}
        for { set j 0 } { $j < [llength $rowdata] } { incr j } {
            set kcal_data [expr { [lindex $rowdata $j]*0.5*1.041308e-21*6.02214e23*0.943*0.943 }] 
            lappend rowdata_kcal $kcal_data
        }
        lappend inthessian_kcal $rowdata_kcal
    }

    return $inthessian_kcal

}
#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::resetDefaultsGenBonded {} {
    # resets the QM settings to the default values
    set ::ForceFieldToolKit::GenBonded::qmProc 1
    set ::ForceFieldToolKit::GenBonded::qmMem 1
    set ::ForceFieldToolKit::GenBonded::qmRoute "MP2"

    # Reset name of output QM file.
    set ::ForceFieldToolKit::GenBonded::com "hess.py"

}
#===========================================================================================================


#===========================================================================================================
# DIHEDRALS PARAMETRIZATION
proc ::ForceFieldToolKit::Psi4::buildFiles_GenDihScan { dihData outPath basename qmProc qmCharge qmMem qmMult qmRoute psf pdb } {

    # assign Psi4 atom names and gather x,y,z for output com file
    ::ForceFieldToolKit::SharedFcns::LonePair::loadMolExcludeLP $psf $pdb
    #mol new $psf; mol addfile $pdb

    # NEW 02/01/2019:
    ::ForceFieldToolKit::SharedFcns::checkElementPDB

    set Gnames {}
    set atom_info {}
    for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
        set temp [atomselect top "index $i"]
        lappend atom_info [list [$temp get element][expr $i+1] [$temp get x] [$temp get y] [$temp get z]]
        lappend Gnames [$temp get element][expr $i+1]
        $temp delete
    }

    # cycle through each dihedral to scan
    foreach sign {1 -1} {
        set scanCount 1
        foreach dih $dihData {
            # change 0-based indices to 1-based
            set zeroInds [lindex $dih 0]
            set oneInds {}
            foreach ind $zeroInds {
                lappend oneInds [expr {$ind + 1}]
            }

            # negative/positive scan
            # open the output file
            if {$sign == 1} {
            set outfile [open ${outPath}/${basename}.scan${scanCount}.pos.py w]
            } elseif {$sign == -1} {
            set outfile [open ${outPath}/${basename}.scan${scanCount}.neg.py w]
            }

            # write the header
            puts $outfile "import psi4"
            puts $outfile "import qcelemental as qcel"
            puts $outfile "import optking"
            puts $outfile "psi4.set_memory(\'$qmMem GB\')"
            puts $outfile "psi4.set_num_threads($qmProc)"
            if {$sign == 1} {
            puts $outfile "psi4.set_output_file('${basename}.scan${scanCount}.pos.out', False)"
            } elseif {$sign == -1} {
            puts $outfile "psi4.set_output_file('${outPath}/${basename}.scan${scanCount}.neg.out', False)"
            }
            puts $outfile ""

            # write coords
            puts $outfile {mol = psi4.geometry("""}
            puts $outfile "$qmCharge $qmMult"
            foreach atom_entry $atom_info {
                puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
            }
            puts $outfile {""")}

            puts $outfile ""
            puts $outfile "xyzs = mol.geometry().np"
            puts $outfile "dihedral = qcel.util.measure_coordinates(xyzs, \[[lindex $zeroInds 0], [lindex $zeroInds 1], [lindex $zeroInds 2], [lindex $zeroInds 3]\], True)"

            puts $outfile ""
            puts $outfile "scan = \[\]"
            puts $outfile "coordinate = \[\]"

            puts $outfile ""
            set stepsize [lindex $dih 2]
            set step [expr int([expr [lindex $dih 1]/$stepsize/2]) + 1]
            puts $outfile "for i in range(0, $step):"
            puts $outfile "    fixD = {\"ranged_dihedral\": \"($oneInds \" + str(dihedral) + \" \" + str(dihedral) + \")\"}"
            puts $outfile "    options = {"
            puts $outfile {        'basis': '6-31g*',}
            puts $outfile {        'mp2_type': 'df',}
            puts $outfile {        'geom_maxiter': 100,}
            puts $outfile {        'dynamic_level': 1,}
            puts $outfile "        }"
            puts $outfile "    psi4.set_options(options)"
            puts $outfile {    json_output = optking.optimize_psi4("mp2", **fixD)}
            puts $outfile {    E = json_output["energies"][-1]}
            puts $outfile "    scan.append((dihedral, E))"
            puts $outfile {    coordinate.append(json_output['trajectory'][-1]['molecule']['geometry'])}
            puts $outfile "    dihedral += [expr $stepsize*$sign]"
            puts $outfile ""
            puts $outfile "print(scan)"

            puts $outfile ""
            puts $outfile {nstep = len(scan)}
            puts $outfile "indices = '[lindex $zeroInds 0] [lindex $zeroInds 1] [lindex $zeroInds 2] [lindex $zeroInds 3]'"
            puts $outfile {energy = json_output["energies"]}

            puts $outfile ""
            if {$sign == 1} {
              set filename_sup $basename.scan${scanCount}.pos.supplement
            } elseif {$sign == -1} {
              set filename_sup $basename.scan${scanCount}.neg.supplement
            }

            puts $outfile {# use indices to write the indices, scan to write energies and dihedrals, and trajectory to write the coordinates}
            puts $outfile "with open('${filename_sup}.out', 'w') as f:"
            puts $outfile {Psi4 dihedral output supplement}
            puts $outfile ""
            puts $outfile {    # write the indices}
            puts $outfile {    f.write("indices \n")}
            puts $outfile {    for i in range(nstep):}
            puts $outfile {        f.write(indices + "\n")}
            puts $outfile {    f.write("end indices \n")}

            puts $outfile {    # write the dihedrals}
            puts $outfile {    f.write("dihedrals \n")}
            puts $outfile {    for i in range(nstep):}
            puts $outfile {        f.write(str(scan[i][0]) + "\n")}
            puts $outfile {    f.write("end dihedrals \n")}

            puts $outfile {    # write the energies}
            puts $outfile {    f.write("energies \n")}
            puts $outfile {    for i in range(nstep):}
            puts $outfile {        f.write(str(scan[i][1]) + "\n")}
            puts $outfile {    f.write("end energies \n")}

            puts $outfile {    # write the coordinates}
            puts $outfile {    f.write("coordinates \n")}
            puts $outfile {    for i in range(nstep):}
            puts $outfile {        s = ""}
            puts $outfile {        for j in range(len(coordinate[i])//3):}
            puts $outfile "            s += \"{\""
            puts $outfile {            for k in range(3):}
            puts $outfile {                s += str(coordinate[i][3*j + k]) + " "}
            puts $outfile "            s += \"} \""
            puts $outfile {        f.write(s + "\n")}
            puts $outfile {    f.write("end coordinates \n")}


        close $outfile

        incr scanCount

        }
    }

}

#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::readLog_TorExplor { log } {

    # initialize some variables
    set indDef {}; set stepEn {}; set stepCoords {}

    # open the log file for reading
    set infile [open $log r]
    while { ![eof $infile] } {
        # read a line at a time
        set inline [string trim [gets $infile]]

        switch -regexp $inline {
            {Initial Parameters} {
                # keep reading until finding the dihedral being scanned
                # and parse out 1-based indices, convert to 0-based
                while { ![regexp {^\!.*D\(([0-9]+),([0-9]+),([0-9]+),([0-9]+)\).*Scan[ \t]+\!$} [string trim [gets $infile]] full_match ind1 ind2 ind3 ind4] && ![eof $infile] } { continue }
                if { [eof $infile] } { return }
                foreach ele [list $ind1 $ind2 $ind3 $ind4] { lappend indDef [expr {$ele - 1}] }
            }

            {Input orientation:} {
                # clear any existing coordinates
                set currCoords {}
                # burn the header
                for {set i 0} {$i<=3} {incr i} { gets $infile }
                # parse coordinates
                while { [string range [string trimleft [set line [gets $infile]] ] 0 0] ne "-" } { lappend currCoords [lrange $line 3 5] }
            }

            {SCF[ \t]*Done:} {
                # parse E(RHF) energy; convert hartrees to kcal/mol
                set currEnergy [expr {[lindex $inline 4] * 627.5095}]
                # NOTE: this value will be overridden if E(MP2) is also found
            }

            {E2.*EUMP2} {
                # convert from Psi4 notation in hartrees to scientific notation
                set currEnergy [expr {[join [split [lindex [string trim $inline] end] D] E] * 627.5095}]
                # NOTE: this overrides the E(RHF) parse from above
            }

            {Optimization completed} {
                # we've reached the optimized conformation
                lappend stepEn $currEnergy
                lappend stepCoords $currCoords
            }

            default {continue}
        }
    }

    close $infile

    # reverse data if it's a negative scan
    if { [regexp {\.neg\.} $log] } {
        set stepEn [lreverse $stepEn]
        set stepCoords [lreverse $stepCoords]
    }

    # return the parsed data
    return [list $indDef $stepEn $stepCoords]
}

#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::parseGlog_DihOpt { debug debugLog GlogFile } {

    # initialize log-wide variables
    set scanDihInds {}; set currDihVal {}; set currCoords {}; set currEnergy {}
    set tmpGlogData {}; set GlogData {}
    set infile [open $GlogFile r]

    # read through the Psi4 output File (Glog)
    while {[eof $infile] != 1} {
        # read a line
        set inline [gets $infile]
        # parse line
        switch -regexp $inline {
            {indices} {
                while { ![string match [string trim [set inline [gets $infile]]] "end indices"] } {
                    lappend scanDihInds $inline
                }
            }
            {dihedral} {
                while { ![string match [string trim [set inline [gets $infile]]] "end dihedrals"] } {
                    lappend currDihVal $inline
                }
            }
            {energies} {
                while { ![string match [string trim [set inline [gets $infile]]] "end energies"] } {
                    lappend currEnergy $inline
                }
            }
            puts $currEnergy
            {coordinates} {
                while { ![string match [string trim [set inline [gets $infile]]] "end coordinates"] } {
                    lappend currCoords $inline
                }
            }
        }; # end of line parse (switch)
    }; # end of cycling through Glog lines (while)

    # for {set i 0} {$i < [llength $scanDihInds]} {incr i}
    for {set i 0} {$i < 3} {incr i} {
        set ltmp [list [lindex $scanDihInds $i] [lindex $currDihVal $i] [lindex $currEnergy $i] [lindex $currCoords $i]]
        lappend tmpGlogData $ltmp
    }

    # if the Psi4 output file runs the scan in negative direction, reverse the order of entries
    # if not explicitely in the negative direction, preserve the order of entries
    if { [lsearch -exact [split $GlogFile \.] "neg"] != -1 } {
        foreach ele [lreverse $tmpGlogData] { lappend GlogData $ele }
    } else {
        foreach ele $tmpGlogData { lappend GlogData $ele }
    }

    # clean up
    close $infile
    return $GlogData

}

#===========================================================================================================
proc ::ForceFieldToolKit::Psi4::resetDefaultsGenDihScan {} {
    # reset QM settings for generation of dihedral scan to the default values
    #
    # set variables
    set ::ForceFieldToolKit::GenDihScan::qmProc 1
    set ::ForceFieldToolKit::GenDihScan::qmCharge 0
    set ::ForceFieldToolKit::GenDihScan::qmMem 1
    set ::ForceFieldToolKit::GenDihScan::qmMult 1
    set ::ForceFieldToolKit::GenDihScan::qmRoute {mp2 6-31g*}

}
#===========================================================================================================
#===========================================================================================================
