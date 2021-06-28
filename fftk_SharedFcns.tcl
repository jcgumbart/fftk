#
# $Id: fftk_SharedFcns.tcl,v 1.18 2020/09/01 16:55:14 johns Exp $
#

namespace eval ::ForceFieldToolKit::SharedFcns {

}
#======================================================
# DEPRECIATED (now that long atom types are supported)
#======================================================
proc ::ForceFieldToolKit::SharedFcns::reTypeFromPSF { psfFile molID } {
    # parses type from PSF and resets appropriately in VMD
    # necessary especially for long atom types (psfgen chokes)

    # short circuit the proc due to depreciation
    puts "reTypeFromPSF has been depreciated, returning without any action"; flush stdout
    return

    set reType {}
    set inFile [open $psfFile r]
    
    # read until !NATOMS section
    while { [lindex [gets $inFile] end] ne "\!NATOM" } {
        continue
    }
    
    # once in NATOMS, read until blank line (end of NATOMS)
    while { [set inline [gets $inFile]] ne "" } {
        lappend reType [lindex $inline 5]
    }
    
    close $inFile

    # reset type
    for {set i 0} {$i < [llength $reType]} {incr i} {
        set temp [atomselect $molID "index $i"]
        $temp set type [lindex $reType $i]
        $temp delete
    }
}
#======================================================
# DEPRECIATED (now that long atom types are supported)
#======================================================
proc ::ForceFieldToolKit::SharedFcns::reChargeFromPSF { psfFile molID } {
    # parses charge from PSF and resets appropriately in VMD
    # necessary especially for long atom types (psfgen chokes) 

    # short circuit the proc due to depreciation
    puts "reTypeFromPSF has been depreciated, returning without any action"; flush stdout
    return

    set reCharge {}
    set inFile [open $psfFile r]
    
    # read until !NATOMS section
    while { [lindex [gets $inFile] end] ne "\!NATOM" } {
        continue
    }
    
    # once in NATOMS, read until blank line (end of NATOMS)
    while { [set inline [gets $inFile]] ne "" } {
        lappend reCharge [lindex $inline 6]
    }
    
    close $inFile

    # reset type and charge 
    for {set i 0} {$i < [llength $reCharge]} {incr i} {
        set temp [atomselect $molID "index $i"]
        $temp set charge [lindex $reCharge $i]
        $temp delete
    }
}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::writeMinConf { name psf pdb parlist {extrabFile ""} } {

   set conf [open "$name.conf" w]
   puts $conf "structure          $psf"    
   puts $conf "coordinates        $pdb"
   puts $conf "paraTypeCharmm     on"
   foreach par $parlist {
      puts $conf "parameters       $par"
   }
   puts $conf "temperature         310"
   puts $conf "exclude             scaled1-4"
   puts $conf "1-4scaling          1.0"
   puts $conf "cutoff              1000.0"
   puts $conf "switching           on"
   puts $conf "switchdist          1000.0"
   puts $conf "pairlistdist        1000.0"
   puts $conf "timestep            1.0 "
   puts $conf "nonbondedFreq       2"
   puts $conf "fullElectFrequency  4 "
   puts $conf "stepspercycle       20"
   puts $conf "outputName          $name"
   puts $conf "restartfreq         1000"
   if { $extrabFile != "" } {
      puts $conf "extraBonds          yes"
      puts $conf "extraBondsFile $extrabFile"
   }
   puts $conf "minimize            1000"
   puts $conf "reinitvels          310"
   puts $conf "run 0"
   close $conf

}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::readParFile { parFile } {
    # reads in a CHARMM parameter file and returns
    # the data as a list
    
    # initialize lists
    set bonds {}
    set angles {}
    set dihedrals {}
    set impropers {}
    set cmaps {}
    set vdws {}
    
    # initialize read state
    set readstate 0
    
    # open the input parameter file
    set inFile [open $parFile r]
    
    # read through it a line at a time
    while { ![eof $inFile] } {
        set inLine [gets $inFile]
        switch -regexp $inLine {
            {^[ \t]*$} { continue }
            {^[ \t]*\*.*} { continue }
            {^[ \t]*!.*} { continue }
            {^[a-z]+} { continue }
            {^BONDS.*} { set readstate BOND }
            {^ANGLES.*} { set readstate ANGLE }
            {^DIHEDRALS.*} { set readstate DIH }
            {^IMPROPER.*} { set readstate IMPROP }
            {^CMAP.*} { set readstate CMAP }
            {^NONBONDED.*} { set readstate VDW }
            {^NBFIX.*} { set readstate 0 }
            {^HBOND.*} { set readstate 0 }
            {^END.*} { break }
            default {
                set prmData [lindex [split $inLine \!] 0]
                set prmComment [string trim [lindex [split $inLine \!] 1]]

                switch -exact $readstate {
                    0 { continue }
                    BOND {
                        #                           { type def }          { k b0 }      { comment }
                        lappend bonds [list [lrange $prmData 0 1] [lrange $prmData 2 3] $prmComment]
                    }
                    ANGLE {
                        #                            { type def }          { k theta }           { kub s }     { comment }
                        lappend angles [list [lrange $prmData 0 2] [lrange $prmData 3 4] [lrange $prmData 5 6] $prmComment]
                    }
                    DIH {
                        #                               { type def }         { k n delta }  { comment }
                        lappend dihedrals [list [lrange $prmData 0 3] [lrange $prmData 4 6] $prmComment]
                    }
                    IMPROP {
                        #                               { type def }                { kpsi                psi0 }     { comment }
                        lappend impropers [list [lrange $prmData 0 3] [list [lindex $prmData 4] [lindex $prmData 6]] $prmComment]
                    }
                    CMAP { continue }
                    VDW {
                        #                         { type def }   { epsilon Rmin/2 }  { eps,1-4 Rmin/2,1-4} { comment }
                        lappend vdws [list [lindex $prmData 0] [lrange $prmData 2 3] [lrange $prmData 5 6] $prmComment]
                    }
                }
                unset prmData; unset prmComment
            }
        }; # end of outer switch
    }; # end while (reading file)
    
    # clean up
    close $inFile
    
    # return
    return [list $bonds $angles $dihedrals $impropers $vdws]

}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::writeParFile { pars filename } {
    # takes in a set of parameters (readcharmmpar format) and a filename
    # writes a charmm-styled parameter file

    # open the output file    
    set outFile [open $filename w]
    
    # HEADER
    # print the header
    #puts $outFile "!*>>>>> CHARMM36 All-Hydrogen Parameter File for Proteins <<<<<<<<"
    #puts $outFile "!*>>>>>>>>>>>>>>>>>>>> and Nucleic Acids <<<<<<<<<<<<<<<<<<<<<<<<<"
    #puts $outFile "!* CHARMM-specific comments to ADM jr. via the CHARMM web site:"
    #puts $outFile "!*                       www.charmm.org"
    #puts $outFile "!*               parameter set discussion forum"
    #puts $outFile "!*"
    #puts $outFile "!"    
    puts $outFile "!============================================================="
    puts $outFile "!"
    puts $outFile "! Parameter file generated by the Force Field ToolKit (ffTK)"
    puts $outFile "!"
    #puts $outFile "! The Force Field ToolKit is a modular set of tools for the"
    #puts $outFile "! development of CHARMM-compatible parameters.  It is"
    #puts $outFile "! available free of charge as a VMD plugin."
    puts $outFile "! For additional information, see:"
    puts $outFile "! http://www.ks.uiuc.edu/Research/vmd/plugins/fftk"
    puts $outFile "! http://www.ks.uiuc.edu/Research/fftk"
    puts $outFile "!"
    puts $outFile "! Authors:"
    puts $outFile "! Christopher G. Mayne"
    puts $outFile "! Beckman Institute for Advanced Science and Technology"
    puts $outFile "! University of Illinois, Urbana-Champaign"
    puts $outFile "! http://www.ks.uiuc.edu/~mayne"
    puts $outFile "! mayne@ks.uiuc.edu"
    puts $outFile "!"
    puts $outFile "! James C. Gumbart"
    puts $outFile "! Georgia Institute of Technology"
    puts $outFile "! http://simbac.gatech.edu"
    puts $outFile "! gumbart_physics.gatech.edu"
    puts $outFile "!"
    #puts $outFile "! Citing the Force Field ToolKit (ffTK)"
    puts $outFile "! If you use parameters developed using ffTK, please cite:"
    #puts $outFile "! C.G. Mayne, J. Saam, K. Schulten, E. Tajkhorshid, J.C. Gumbart. J. Comput. Chem. 2013, DOI: 10.1002/jcc.23422."
    puts $outFile "! C.G. Mayne, J. Saam, K. Schulten, E. Tajkhorshid, J.C. Gumbart. J. Comput. Chem. 2013, 34, 2757-2770."
    puts $outFile "!"
    puts $outFile "!============================================================="
    
    # BONDS
    # determine the maximum field widths to make it look pretty
    set b1max 0
    set b2max 0
    foreach bondDef [lindex $pars 0] {
        set b1l [string length [lindex $bondDef 0 0]]
        set b2l [string length [lindex $bondDef 0 1]]
        if {$b1l > $b1max} {
            set b1max $b1l
        }
        if {$b2l > $b2max} {
            set b2max $b2l
        }
    }
    
    # print the bonds section
    puts $outFile "\nBONDS"
    puts $outFile "!V(bond) = Kb(b - b0)**2"
    puts $outFile "!"
    puts $outFile "!Kb: kcal/mole/A**2"
    puts $outFile "!b0: A"
    puts $outFile "!"
    puts $outFile "!atom type Kb b0"
    puts $outFile "!"
    foreach bondDef [lindex $pars 0] {
        set at1 [lindex $bondDef 0 0]
        set at2 [lindex $bondDef 0 1]
        set k [lindex $bondDef 1 0]
        set b [lindex $bondDef 1 1]
        set comment [lindex $bondDef 2]
        puts $outFile "[format %-*s $b1max $at1]  [format %-*s $b2max $at2]  [format %-9s [format %.3f $k]]  [format %-7s [format %.3f $b]]  ! $comment"
    }
    
    # ANGLES
    # determine the maximum field widths to make it look pretty
    set a1max 0
    set a2max 0
    set a3max 0
    foreach angleDef [lindex $pars 1] {
        set a1l [string length [lindex $angleDef 0 0]]
        set a2l [string length [lindex $angleDef 0 1]]
        set a3l [string length [lindex $angleDef 0 2]]
        if {$a1l > $a1max} {
            set a1max $a1l
        }
        if {$a2l > $a2max} {
            set a2max $a2l
        }
        if {$a3l > $a3max} {
            set a3max $a3l
        }
    }
    
    # print the angles section
    puts $outFile "\nANGLES"
    puts $outFile "!"
    puts $outFile "!V(angle) = Ktheta(Theta - Theta0)**2"
    puts $outFile "!"
    puts $outFile "!V(Urey-Bradley) = Kub(S - S0)**2"
    puts $outFile "!"
    puts $outFile "!Ktheta: kcal/mole/rad**2"
    puts $outFile "!Theta0: degrees"
    puts $outFile "!Kub: kcal/mole/A**2 (Urey-Bradley)"
    puts $outFile "!S0: A"
    puts $outFile "!"
    puts $outFile "!atom types     Ktheta    Theta0   Kub     S0"
    puts $outFile "!"
    puts $outFile "!"
    foreach angleDef [lindex $pars 1] {
        set at1 [lindex $angleDef 0 0]
        set at2 [lindex $angleDef 0 1]
        set at3 [lindex $angleDef 0 2]
        set ktheta [lindex $angleDef 1 0]
        set theta [lindex $angleDef 1 1]
        set kub [lindex $angleDef 2 0]
        if { $kub ne ""} { set kub "[format %-7s [format %.2f $kub]]" }
        set s [lindex $angleDef 2 1]
        if { $s ne "" } { set s "[format %-7s [format %.4f $s]]" }
        set comment [lindex $angleDef 3]
        puts $outFile "[format %-*s $a1max $at1]  [format %-*s $a2max $at2]  [format %-*s $a3max $at3]  [format %-7s [format %.3f $ktheta]]  [format %-7s [format %.3f $theta]]  $kub  $s ! $comment"
    }
    
    # DIHEDRALS
    # determine the maximum field widths to make it look pretty
    set d1max 0
    set d2max 0
    set d3max 0
    set d4max 0
    foreach dihDef [lindex $pars 2] {
        set d1l [string length [lindex $dihDef 0 0]]
        set d2l [string length [lindex $dihDef 0 1]]
        set d3l [string length [lindex $dihDef 0 2]]
        set d4l [string length [lindex $dihDef 0 3]]
        if {$d1l > $d1max} {
            set d1max $d1l
        }
        if {$d2l > $d2max} {
            set d2max $d2l
        }
        if {$d3l > $d3max} {
            set d3max $d3l
        }
        if {$d4l > $d4max} {
            set d4max $d4l
        }
    }  
    
    # print the dihedrals section
    puts $outFile "\nDIHEDRALS"
    puts $outFile "!"
    puts $outFile "!V(dihedral) = Kchi(1 + cos(n(chi) - delta))"
    puts $outFile "!"
    puts $outFile "!Kchi: kcal/mole"
    puts $outFile "!n: multiplicity"
    puts $outFile "!delta: degrees"
    puts $outFile "!"
    puts $outFile "!atom types             Kchi    n   delta"
    puts $outFile "!"
    foreach dihDef [lindex $pars 2] {
        set at1 [lindex $dihDef 0 0]
        set at2 [lindex $dihDef 0 1]
        set at3 [lindex $dihDef 0 2]
        set at4 [lindex $dihDef 0 3]
        set k [lindex $dihDef 1 0]
        set n [lindex $dihDef 1 1]
        set delta [lindex $dihDef 1 2]
        set comment [lindex $dihDef 2]
        puts $outFile "[format %-*s $d1max $at1]  [format %-*s $d2max $at2]  [format %-*s $d3max $at3]  [format %-*s $d4max $at4]  [format %-7s [format %.4f $k]]  [format %-1s $n]  [format %-5s [format %.2f $delta]] ! $comment"
    }
    
    # IMPROPERS
    # determine the maximum field widths to make it look pretty
    set i1max 0
    set i2max 0
    set i3max 0
    set i4max 0
    foreach imprDef [lindex $pars 3] {
        set i1l [string length [lindex $imprDef 0 0]]
        set i2l [string length [lindex $imprDef 0 1]]
        set i3l [string length [lindex $imprDef 0 2]]
        set i4l [string length [lindex $imprDef 0 3]]
        if {$i1l > $i1max} {
            set i1max $i1l
        }
        if {$i2l > $i2max} {
            set i2max $i2l
        }
        if {$i3l > $i3max} {
            set i3max $i3l
        }
        if {$i4l > $i4max} {
            set i4max $i4l
        }        
    }
    
    # print the impropers section
    puts $outFile "\nIMPROPER"
    puts $outFile "!"
    puts $outFile "!V(improper) = Kpsi(psi - psi0)**2"
    puts $outFile "!"
    puts $outFile "!Kpsi: kcal/mole/rad**2"
    puts $outFile "!psi0: degrees"
    puts $outFile "!note that the second column of numbers (0) is ignored"
    puts $outFile "!"
    puts $outFile "!atom types           Kpsi                   psi0"
    puts $outFile "!"
    foreach imprDef [lindex $pars 3] {
        set at1 [lindex $imprDef 0 0]
        set at2 [lindex $imprDef 0 1]
        set at3 [lindex $imprDef 0 2]
        set at4 [lindex $imprDef 0 3]
        set kpsi [lindex $imprDef 1 0]
        set n 0
        set psi0 [lindex $imprDef 1 1]
        if { $psi0 ne "" } { set psi0 [format %.2f $psi0] }
        set comment [lindex $imprDef 2]
        puts $outFile "[format %-*s $d1max $at1]  [format %-*s $d2max $at2]  [format %-*s $d3max $at3]  [format %-*s $d4max $at4]  [format %-7s [format %.4f $kpsi]]  [format %-1s $n]  [format %-5s $psi0]  ! $comment"
    }
    
    # NONBONDED
    # determine maximum field widths to make it look pretty
    set atmax 0
    foreach nonbDef [lindex $pars 4] {
        set atl [string length [lindex $nonbDef 0]]
        if {$atl > $atmax } {
            set atmax $atl
        }
    }
    
    # print the nonbonded section
    puts $outFile "\nNONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -"
    puts $outFile "cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5 "
    puts $outFile "!"
    puts $outFile "!V(Lennard-Jones) = Eps,i,j\[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6\]"
    puts $outFile "!"
    puts $outFile "!epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)"
    puts $outFile "!Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j"
    puts $outFile "!"
    puts $outFile "!atom  ignored    epsilon      Rmin/2   ignored   eps,1-4       Rmin/2,1-4"
    puts $outFile "!"
    foreach nonbDef [lindex $pars 4] {
        # for whatever reason, ::Pararead::read_charmm_parameters appends the
        # END statement to non-bonded parameters, which throws an error
        if { [lindex $nonbDef 0] eq "END" } { continue }
        set at [lindex $nonbDef 0]
        set eps [lindex $nonbDef 1 0]
        set rmin [lindex $nonbDef 1 1]
        set eps2 [lindex $nonbDef 2 0]
        if { $eps2 ne "" } { set eps2 "0.0  [format %-8s [format %.6f $eps2]]" }
        set rmin2 [lindex $nonbDef 2 1]
        if { $rmin2 ne "" } { set rmin2 "[format %-8s [format %.6f $rmin2]]" }
        set comment [lindex $nonbDef 3]
        puts $outFile "[format %-*s $atmax $at]  0.0  [format %-8s [format %.6f $eps]]  [format %-8s [format %.6f $rmin]]  $eps2  $rmin2  ! $comment"
    }
    
    # WRAP UP
    puts $outFile "\nEND"
    
    close $outFile    
}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::avgZmatReplicates { molID zmat } {
    # averages replicate zmatrix entries
    
    # passed: molid for a properly (re)typed molecule, zmatrix
    # returns: reformatted version of zmatrix with averaged replicate values (replicates are defined by identical typedef)
    
    # initialize arrays
    array unset bonds; array set bonds {}; # entry format: bonds(typeDef) = { {cumulative k} {cumulative b0} {count} }
    array unset angles; array set angles {}; # entry format: angles(typeDef) = { {cumulative k} {cumulative theta} {count} }
    array unset dihedrals; array set dihedrals {}; # entry format: dihedrals(typeDef) = { {cumulative k} {last n} {cumulative phase shift} {count} }
    
    # burn the zmat header
    set zmat [lreplace $zmat 0 0]
    
    # process each zmat entry
    foreach zmatEntry $zmat {

        # common processing    
        # convert indDef to typeDef
        set typeDef {}
        foreach ind [lindex $zmatEntry 2] {
           lappend typeDef [[atomselect $molID "index $ind"] get type]
        }

        # bond/angle/dihedral-specific processing
        switch -exact [lindex $zmatEntry 1] {
            {bond} {
                # check to see if fwd or rev typeDef has already been added to the array
                set testfwd [info exists bonds($typeDef)]
                set testrev [info exists bonds([lreverse $typeDef])]
                if { $testfwd == 0 && $testrev == 0 } {
                    # this is a new typeDef; build an entry
                    set bonds($typeDef) [list [lindex $zmatEntry 4 0] [lindex $zmatEntry 4 1] 1]
                } else {
                    # entry is a typeDef replicate; update accordingly                    
                    if { $testfwd } { set matchedDef $typeDef } else { set matchedDef [lreverse $typeDef] }
                    lset bonds($matchedDef) 0 [expr { [lindex $bonds($matchedDef) 0] + [lindex $zmatEntry 4 0]}]
                    lset bonds($matchedDef) 1 [expr { [lindex $bonds($matchedDef) 1] + [lindex $zmatEntry 4 1]}]
                    lset bonds($matchedDef) 2 [expr { [lindex $bonds($matchedDef) 2] + 1}]
                }
            }
            {lbend} -
            {angle} {
                # check to see if fwd or rev typeDef has already been added to the array
                set testfwd [info exists angles($typeDef)]
                set testrev [info exists angles([lreverse $typeDef])]
                if { $testfwd == 0 && $testrev == 0 } {
                    # this is a new typeDef; build an entry
                    set angles($typeDef) [list [lindex $zmatEntry 4 0] [lindex $zmatEntry 4 1] 1]
                } else {
                    # entry is a typeDef replicate; update accordingly                    
                    if { $testfwd } { set matchedDef $typeDef } else { set matchedDef [lreverse $typeDef] }
                    lset angles($matchedDef) 0 [expr { [lindex $angles($matchedDef) 0] + [lindex $zmatEntry 4 0]}]
                    lset angles($matchedDef) 1 [expr { [lindex $angles($matchedDef) 1] + [lindex $zmatEntry 4 1]}]
                    lset angles($matchedDef) 2 [expr { [lindex $angles($matchedDef) 2] + 1}]
                }
            }
            {dihed} {
                # there are problems with zmatqmEff in the bonds/angles opt routine which causes errors here
                # therefore i need to skip this section until that is fixed
                continue
                # check to see if fwd or rev typeDef has already been added to the array
                set testfwd [info exists dihedrals($typeDef)]
                set testrev [info exists dihedrals([lreverse $typeDef])]
                if { $testfwd == 0 && $testrev == 0 } {
                    # this is a new typeDef; build an entry
                    set dihedrals($typeDef) [list [lindex $zmatEntry 4 0] [lindex $zmatEntry 4 1] [lindex $zmatEntry 4 2] 1]
                } else {
                    # entry is a typeDef replicate; update accordingly                    
                    if { $testfwd } { set matchedDef $typeDef } else { set matchedDef [lreverse $typeDef] }
                    lset dihedrals($matchedDef) 0 [expr { [lindex $dihedrals($matchedDef) 0] + [lindex $zmatEntry 4 0]}]
                    lset dihedrals($matchedDef) 1 [lindex $zmatEntry 4 1]
                    lset dihedrals($matchedDef) 2 [expr { [lindex $dihedrals($matchedDef) 2] + [lindex $zmatEntry 4 2]}]
                    lset dihedrals($matchedDef) 3 [expr { [lindex $dihedrals($matchedDef) 3] + 1}]
                }
            }
            default { continue }
        }; # end switch
    }; # end zmat element processing
    
    
    # rebuild a pseudo-zmat (reformatted)
    set returnZmat {}
    # bonds
    foreach typeDef [array names bonds] {
        set avgK [format %.4f [expr { [lindex $bonds($typeDef) 0] / [lindex $bonds($typeDef) 2] }]]
        set avgB0 [format %.4f [expr { [lindex $bonds($typeDef) 1] / [lindex $bonds($typeDef) 2] }]]
        lappend returnZmat [list bond $typeDef $avgK $avgB0]
    }
    # angles
    foreach typeDef [array names angles] {
        set avgK [format %.4f [expr { [lindex $angles($typeDef) 0] / [lindex $angles($typeDef) 2] }]]
        set avgTheta [format %.4f [expr { [lindex $angles($typeDef) 1] / [lindex $angles($typeDef) 2] }]]
        lappend returnZmat [list angle $typeDef $avgK $avgTheta]
    }
    # dihedrals
    foreach typeDef [array names dihedrals] {
        set avgK [format %.4f [expr { [lindex $dihedrals($typeDef) 0] / [lindex $dihedrals($typeDef) 3] }]]
        set n [lindex $dihedrals($typeDef) 1]
        set avgDelta [format %.4f [expr { [lindex $dihedrals($typeDef) 2] / [lindex $dihedrals($typeDef) 3] }]]
        lappend returnZmat [list dihed $typeDef $avgK $n $avgDelta]
    }    
    
    # return the reformatted pseudo-zmat
    return $returnZmat    
}
#======================================================
#======================================================
# ParView
#======================================================
namespace eval ::ForceFieldToolKit::SharedFcns::ParView {
    variable objList
    array unset objList; array set objList {}
}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::ParView::clearParViewObjList {args} {
    # clears the objects and objList for a specified molid
    # localize relevant variables
    variable objList

    # args defaults
    set molid [molinfo top]

    # parse passed args
    foreach {flag val} $args { set [string range $flag 1 end] $val }

    # delete any existing objects and wipe the objList
    if { [info exists objList($molid)] } {
        foreach obj $objList($molid) { graphics $molid delete $obj }        
    }
    set objList($molid) {}
}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::ParView::addColorObj {args} {
    # adds a color object to the objList
    # localize relevant variables
    variable objList

    # args defaults
    set molid [molinfo top]
    set color blue

    # parse passed args
    foreach {flag val} $args { set [string range $flag 1 end] $val }

    # set the color and add it to the objList element
    lappend objList($molid) [graphics $molid color $color]
}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::ParView::addMaterialObj {args} {
    # adds a color object to the objList
    # localize relevant variables
    variable objList

    # args defaults
    set molid [molinfo top]
    set material Diffuse

    # parse passed args
    foreach {flag val} $args { set [string range $flag 1 end] $val }

    # set the color and add it to the objList element
    lappend objList($molid) [graphics $molid material $material]
}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::ParView::addParObject {args} {
    # creates a graphic object for a given type and adds it to the objList
    # localize relevant variables
    variable objList
    
    # argument flags
    # -indList = list of inds
    # -frame = now|end|<integer>
    # -type = atom|bond|angle|dihedral|improper
    # -molid = molecule id

    # set default arguments
    set frame 0
    set indices {}
    set type {}
    set molid [molinfo top]
    set radius 0.1
    set resolution 30

    # parse passed args
    foreach {flag val} $args { set [string range $flag 1 end] $val }

    # construct a list of coordinates from the index list
    set xyzList {}
    foreach ind $indices {
        set sel [atomselect $molid "index $ind" frame $frame]
        lappend xyzList [lindex [$sel get {x y z}] 0]
        $sel delete
    }

    # scale elements to create layers that render properly with tachyon AO
    set bondRadius $radius
    set angRadius [expr {$radius+0.001}]
    set dihRadius [expr {$radius+0.002}]
    set imprpRadius [expr {$radius}]

    # based on type value, take action
    switch $type {
        {atom} {
            # construct a sphere
            lappend objList($molid) [graphics $molid sphere [lindex $xyzList 0] radius [expr {2*$radius}] resolution $resolution]
        }

        {bond} {
            # construct a sphere at each coord
            foreach xyz $xyzList { lappend objList($molid) [graphics $molid sphere $xyz radius [expr {1.4*$bondRadius}] resolution $resolution] }
            # construct a cylinder from 0->1
            lappend objList($molid) [graphics $molid cylinder [lindex $xyzList 0] [lindex $xyzList 1] radius $bondRadius resolution $resolution filled yes]
            # construct a middle-cylinder to allow viz of bond when other pars are also shown
            set vMiddle [vecscale 0.5 [vecadd [lindex $xyzList 0] [lindex $xyzList 1]]]
            set vLow [vecadd [vecscale 0.5 [lindex $xyzList 0]] [vecscale 0.5 $vMiddle]]
            set vHigh [vecadd [vecscale 0.5 $vMiddle] [vecscale 0.5 [lindex $xyzList 1]]]
            lappend objList($molid) [graphics $molid cylinder $vLow $vHigh radius [expr {2*$radius}] resolution $resolution filled yes]
        }

        {angle} {
            # construct a sphere at each coord
            foreach xyz $xyzList { lappend objList($molid) [graphics $molid sphere $xyz radius [expr {1.4*$angRadius}] resolution $resolution] }
            # construct cylinders from 0->1->2
            lappend objList($molid) [graphics $molid cylinder [lindex $xyzList 0] [lindex $xyzList 1] radius $angRadius resolution $resolution filled yes] 
            lappend objList($molid) [graphics $molid cylinder [lindex $xyzList 1] [lindex $xyzList 2] radius $angRadius resolution $resolution filled yes]
            # construct a triangle in the angle armpit
            set vMiddle1 [vecadd [vecscale 0.4 [lindex $xyzList 0]] [vecscale 0.6 [lindex $xyzList 1]]]
            set vMiddle2 [vecadd [vecscale 0.6 [lindex $xyzList 1]] [vecscale 0.4 [lindex $xyzList 2]]]
            lappend objList($molid) [graphics $molid triangle $vMiddle1 [lindex $xyzList 1] $vMiddle2]
        }

        {dihedral} {
            # construct a sphere at each coord
            foreach xyz $xyzList { lappend objList($molid) [graphics $molid sphere $xyz radius [expr {1.4*$dihRadius}] resolution $resolution] }
            # construct cylinders from 0->1->2->3
            lappend objList($molid) [graphics $molid cylinder [lindex $xyzList 0] [lindex $xyzList 1] radius $dihRadius resolution $resolution filled yes] 
            lappend objList($molid) [graphics $molid cylinder [lindex $xyzList 1] [lindex $xyzList 2] radius $dihRadius resolution $resolution filled yes]
            lappend objList($molid) [graphics $molid cylinder [lindex $xyzList 2] [lindex $xyzList 3] radius $dihRadius resolution $resolution filled yes] 
            # construct the torsion-denoting cylinder in between 1->2
            set vMiddle [vecscale 0.5 [vecadd [lindex $xyzList 1] [lindex $xyzList 2]]]
            #set vLow [vecscale 0.5 [vecadd [lindex $xyzList 1] $vMiddle]]
            set vLow [vecadd [vecscale 0.1 [lindex $xyzList 1]] [vecscale 0.9 $vMiddle]]
            #set vHigh [vecscale 0.5 [vecadd $vMiddle [lindex $xyzList 2]]]
            set vHigh [vecadd [vecscale 0.9 $vMiddle] [vecscale 0.1 [lindex $xyzList 2]]]
            lappend objList($molid) [graphics $molid cylinder $vLow $vHigh radius [expr {3*$radius}] resolution $resolution filled yes]
        }

        {improper} {
            # construct a sphere at each coord
            foreach xyz $xyzList { lappend objList($molid) [graphics $molid sphere $xyz radius [expr {1.4*$imprpRadius}] resolution $resolution] }
            # cylinder from 0->1, 0->2, 0->3
            lappend objList($molid) [graphics $molid cylinder [lindex $xyzList 0] [lindex $xyzList 1] radius $imprpRadius resolution $resolution filled yes] 
            lappend objList($molid) [graphics $molid cylinder [lindex $xyzList 0] [lindex $xyzList 2] radius $imprpRadius resolution $resolution filled yes] 
            lappend objList($molid) [graphics $molid cylinder [lindex $xyzList 0] [lindex $xyzList 3] radius $imprpRadius resolution $resolution filled yes]
            # construct a short cylinder (disc) around center extending perpendicular to plane of improper
            # FIXME 
        }
    }
}
#======================================================
namespace eval ::ForceFieldToolKit::SharedFcns::LonePair {
    variable LPinfo
}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::LonePair::initFromDict { molID LPDict } {
    # initialize the lp list by using a dictionary provided. Generate the nonexisting LPs
    variable LPinfo

    set LPinfo {}
    set toBuildLP {}

    if {[dict size $LPDict] <= 0} {
        return $molID
    }

    # build dict with index as key if defined
    dict for {name value} $LPDict {
        dict with value {
            if {$index >= 0} {
                dict set LPinfo $index name   $name
                dict set LPinfo $index type   $type
                dict set LPinfo $index charge $charge
                dict set LPinfo $index host1  $host1
                dict set LPinfo $index host2  $host2
                dict set LPinfo $index dist   $dist
            } else {
                dict set toBuildLP $name $value
            }
        }
    }

    # add lone pair atoms if needed
    set numToBuildLP [dict size $toBuildLP]
    if {$numToBuildLP > 0} {
        set lpMol [mol new atoms $numToBuildLP]
        animate dup $lpMol

        set n [molinfo $molID get numatoms]
        set i 0
        dict for {name info} $toBuildLP {
            set newindex [expr $n + $i]
            set sel [atomselect $lpMol "index $i"]
            dict with info {
                dict set LPinfo $newindex name   $name
                dict set LPinfo $newindex type   $type
                dict set LPinfo $newindex charge $charge
                dict set LPinfo $newindex host1  $host1
                dict set LPinfo $newindex host2  $host2
                dict set LPinfo $newindex dist   $dist

                $sel set name   $name
                $sel set type   $type
                $sel set charge $charge
                $sel set mass   0
                if [info exists penalty] {
                    $sel set beta $penalty
                }
            }
            $sel delete
            incr i
        }
        set sel0 [atomselect $molID all]
        set sel1 [atomselect $lpMol all]
        set betas "[$sel0 get beta] [$sel1 get beta]"
        $sel0 delete
        $sel1 delete

        set mergeMol [::TopoTools::mergemols "$molID $lpMol"]
        set sel [atomselect $mergeMol all]
        $sel set beta $betas
        $sel delete

        mol delete $molID
        mol delete $lpMol
        set molID $mergeMol
    }

    # construct lp positions
    set selnolp [atomselect $molID "mass > 0"]
    set selall  [atomselect $molID all]
    $selall set {x y z} [::ForceFieldToolKit::SharedFcns::LonePair::addLPCoordinate [$selnolp get {x y z}]]

    # add bond
    dict for {index value} $LPinfo {
        dict with value {
            topo addbond $index $host1
        }
    }

    return $molID

}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::LonePair::initFromPSF { psf {resNameLimit ""} } {
    # initialize the lp list and its hosts from psf

    variable LPinfo

    set LPinfo {}

    set fp [open $psf r]
    set lines [split [read $fp] "\n"]
    close $fp

    set distances {}
    set nhosts {}
    set nrec -1
    set mode header
    foreach line $lines {
        switch $mode {
            header {
                # find if the psffile is in extended mode
                if {[lindex $line 1] eq "EXT"} {
                    set charmmext true
                } else {
                    set charmmext false
                }
                set mode findatom
            }
            findatom {
                # Look for the NATOM
                if {[string first "NATOM" $line] >= 0} {
                    set nrec [lindex $line 0]
                    set index 0

                    set mode readatom
                }
            }
            readatom {
                # read atom info
                set rc [scan $line "%d %8s %8s %8s %8s %8s %lf %lf %d" \
                    num segname resid resname name atype charge mass lpd]

                # psfgen checked if an atom is lone pair by checking if $lpd is -1, 
                # but NAMD disregard it and check if atom is massless
                # Sticking with NAMD's interpretation

                # Mass: smallest possible mass reresented in psf
                # Resname: equal to the requirement if set
                if {$mass < 0.0001 && ("$resNameLimit" eq "" || "$resname" eq "$resNameLimit")} {
                    set lpindex [expr $num - 1]
                    dict set LPinfo $lpindex name   $name
                    dict set LPinfo $lpindex type   $atype
                    dict set LPinfo $lpindex charge $charge
                }

                incr index
                if {$index >= $nrec} {
                    set mode findnumlp
                }
            }
            findnumlp {
                # Look for the NUMLP
                if {[string first "NUMLP" $line] >= 0} {
                    set nrec [lindex $line 0]
                    set nent [lindex $line 1]
                    set index 0

                    set mode readnumlp
                }
            }
            readnumlp {
                # read lonepair info
                set rc [scan $line "%d %d %s %f %f %f" nhost lpindex type dist ang dih]
                if {$rc != 6} {
                    tk_messageBox -type ok -icon warning -message "PSF file parsing error!" -detail "$line"
                    return
                }

                # print warning if not colinear LP
                if {$nhost != 2} {
                    puts "WARNING: non colinear LP are not supported by ffTK!"
                }

                lappend distances $dist
                lappend nhosts    $nhost

                incr index
                if {$index >= $nrec} {
                    set etype 0
                    set index 0
                    set entry_i 0

                    set mode readLPEntries
                }
            }
            readLPEntries {
                # each line has 8 entries max
                for {set i 0} {$i < 8} {incr i} {
                    set x [expr [lindex $line $i] - 1]
                    switch $etype {
                        0 {
                            if {[dict exists $LPinfo $x]} {
                                set lpindex $x
                                dict set LPinfo $lpindex dist [lindex $distances $index]
                            } else {
                                set lpindex ""
                            }
                        }
                        1 {
                            if {$lpindex ne ""} {
                                dict set LPinfo $lpindex host1 $x
                            }
                        }
                        2 {
                            if {$lpindex ne ""} {
                                dict set LPinfo $lpindex host2 $x
                            }
                        }
                        default {
                            # do nth for host3 or more
                        }
                    }

                    incr entry_i
                    if {$entry_i >= $nent} {
                        set mode quit
                        break
                    }

                    # if read all entries for this lp, reset
                    if {$etype >= [lindex $nhosts $index]} {
                        set etype 0
                    } else {
                        incr etype
                    }
                }
            }
            quit {
                break
            }
        }
    }

    return [mol new $psf]
}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::LonePair::isLP { index } {
    # return if the atom is a lone pair particle

    variable LPinfo
    if {[info exists LPinfo]} {
        return [dict exists $LPinfo $index]
    } else {
        return false
    }
}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::LonePair::numLP { } {
    # return the number of lonepairs

    variable LPinfo
    if {[info exists LPinfo]} {
        return [dict size $LPinfo]
    } else {
        return 0
    }
}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::LonePair::addLPCoordinate { coords } {
    # generate lone pair position on QMcoords

    variable LPinfo

    dict for {index value} $LPinfo {
        dict with value {
            set xyz_h1 [lindex $coords $host1]
            set xyz_h2 [lindex $coords $host2]

            set dir [vecnorm [vecsub $xyz_h1 $xyz_h2]]
            set pos [vecadd $xyz_h1 [vecscale $dist $dir]]

            set coords [linsert $coords $index $pos]
        }
    }

    return $coords
}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::LonePair::loadMolExcludeLP { psf pdb } {
    # generate a psf wihtout lone pair
    
    set molID [mol new $psf]
    mol addfile $pdb
    set all [atomselect $molID all]
    set nolp [atomselect $molID "mass > 0"]

    if {[$all num] > [$nolp num]} {
        $nolp writepsf "nolp.psf"
        $nolp writepdb "nolp.pdb"

        mol delete $molID
        set molID [mol new "nolp.psf"]
        mol addfile "nolp.pdb"

        file delete "nolp.psf"
        file delete "nolp.pdb"
    }

    $all delete
    $nolp delete

    return $molID
}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::LonePair::writePSF { molID fname } {
    variable LPinfo

    set sel [atomselect $molID all]
    $sel writepsf $fname
    $sel delete

  # By Brian Radak (copied from psfgen: topo_mol_output.c)
  # Requirements for NAMD as I understand them:
  # 1) LP has 0 mass and follows parent atoms (the "LP hosts", in this case CL and C6)
  #
  # 2) bond between LP and parent should be optional in the bonds section, but
  #    this might be necessary for building migration groups properly?
  #
  # 3) no angles, dihedrals, etc. should contain LP, this might just be because those
  #    types are not defined

  # 4) the new "!NUMLP NUMLPH" section in the PSF should exist. The presence of
  #    this section automatically toggles the now deprecated "lonepairs on"
  #    keyword and is also required for "drude on". Format is as follows:
  # =========
  # <# of lonepairs> <# of lonepairs + # of lphosts> !NUMLP NUMLPH
  # .
  # :
  # <# of lphosts for this lonepair> <the pointer for this lonepair> F <distance> <angle> <dihedral>
  # .
  # :
  # <lonepair indices>
  # <lonepair index> <lphost1> <lphost2> .... <lonepair index> ...
  # ==========

  # Note the mislabeling of NUMLPH in the comment, the actual number of hosts is
  # the second number minus the first.

  # For each entry, a collinear lonepair has 2 hosts and a bisectory has 3. I
  # forget what F mean ("fixed"?) this is the only option I've seen and the
  # only one NAMD accepts. All entries require a <distance>, <angle>, <dihedral>
  # specification. For collinear lonepairs dihedral is read, but ignored. <angle>
  # is interpreted as a scale parameter that shifts the origin for the colinear x distance
  # back along the bond between the parent atom and the other host. A value of 0.0 means
  # that the origin of the distance is the parent atom.
  #
  # The lonepairs are essentially double indexed (starting at one in Fortran style),
  # once as an atom and again within the list of lonepairs and hosts. So for the example,
  # the first lonepair has index 1, which refers to the 13, which is the lonepair
  # atom index. The 2 for that entry indicates two extra host entries (12 and 11,
  # lonepair indices 2 and 3) and because it is collinear, the lonepair is bound to atom 12.

  # An additional lonepair would have index 4, add either 2 or 3 more hosts to the
  # index list and thus add 3 or 4 more entries. The index list should wrap to the
  # next line after every 8 entries.
  # Example, lets add two 5 point SWM4 waters with bisector lonepairs after the chorobenzene.
  # The new section would be:

  # 3        11 !NUMLP NUMLPH
  # 2         1   F   1.64000       0.00000       0.00000
  # 3         4   F  -0.240345     0.00000       0.00000
  # 3         8   F  -0.240345     0.00000       0.00000
  # 13        12        11        16       14      17        18        21
  # 19        22        23

    if { [info exist LPinfo] && [dict size $LPinfo] > 0 } {
        puts "LPinfo exist! [dict size $LPinfo]"
        set fp_in [open $fname r]
        set lines [split [read -nonewline $fp_in] "\n"]
        close $fp_in

        # get lone pair indices
        set lp_indices [dict keys $LPinfo]
        set lp_indices [lsort -integer $lp_indices]

        set fp_out [open "tmp.psf" w]
        set mode header
        foreach line $lines {
            switch $mode {
                header {
                    # copy line
                    puts $fp_out $line

                    # find if the psffile is in extended mode
                    if {[lindex $line 1] eq "EXT"} {
                        set charmmext true
                    } else {
                        set charmmext false
                    }
                    set mode findatom
                }
                findatom {
                    # copy line
                    puts $fp_out $line

                    # Look for the NATOM
                    if {[string first "NATOM" $line] >= 0} {
                        set nrec [lindex $line 0]
                        set index 0
                        set lpi 0

                        set mode readatom
                    }
                }
                readatom {
                    # read atom info
                    set rc [scan $line "%d %8s %8s %8s %8s %8s %lf %lf %d" \
                        num segname resid resname name atype charge mass lpd]

                    if {$num-1 == [lindex $lp_indices $lpi]} {
                        # if lone pair, change lpd to -1
                        set lpd -1

                        if {$charmmext} {
                            set fmt "%10d %-8s %-8s %-8s %-8s %-6s %10.6f    %10.4f  %10d"
                        } else {
                            set fmt "%8d %-4s %-4s %-4s %-4s %-4s %10.6f    %10.4f  %10d"
                        }
                        puts $fp_out [format $fmt $num $segname $resid $resname $name $atype $charge $mass $lpd]

                        incr lpi
                    } else {
                        # if not lone pair, just copy
                        puts $fp_out $line
                    }

                    incr index
                    if {$index >= $nrec} {
                        set mode findnumlp
                    }
                }
                findnumlp {
                    # Look for the NUMLP
                    if {[string first "NUMLP" $line] >= 0} {
                        set nrec [lindex $line 0]
                        set index 0

                        # discard the existing NUMLP entries
                        # in case in the future atomselect can write NUMLP entries
                        set mode discard
                    } else {
                        # copy line
                        puts $fp_out $line
                    }
                }
                discard {
                    incr index
                    if {$index >= $nrec} {
                        set mode default
                    }
                }
                default {
                    # copy line
                    puts $fp_out $line
                }
            }
        }


        # currently we only handle colinear lps
        # so numlph = numlp*3
        set n [dict size $LPinfo]
        puts $fp_out [format "%8d %8d !NUMLP NUMLPH" $n [expr $n*3]]

        set i 1
        dict for {index info} $LPinfo {
            dict with info {
                # 2 hosts, index, F, distance, angle, dihedral
                puts $fp_out [format "%8d %8d %-6s %10.4f %10.4f %10.4f" 2 $i "F" $dist 0 0]
                incr i 3
            }
        }

        set i 1
        dict for {index info} $LPinfo {
            dict with info {
                foreach var {index host1 host2} {
                    # newline every 8 entries
                    if {$i >= 8} {
                        puts $fp_out ""
                        set i 0
                    }

                    puts -nonewline $fp_out [format " %9d" [eval expr $$var + 1]]
                    incr i
                }
            }
        }
        puts $fp_out ""
        puts $fp_out ""

        close $fp_out

        file rename -force "tmp.psf" $fname
    }
}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::checkWhichQM { outFile } {
    # Check if the QM output file was generated by Gaussian or ORCA
    
    variable qmSoft $::ForceFieldToolKit::qmSoft
    variable NameQM "Not recognized"
    set exitVar 0
    set inFile [open $outFile r]

    # read and match
    while { ![eof $inFile] } { 
	set line [string trim [string map { \" {} } [gets $inFile]]]
	# Check if Gaussian
	if { [string match -nocase "*Gaussian, Inc.*" $line] } {
	       set NameQM "Gaussian"
               break 
	# Check if ORCA
	} elseif { [string match -nocase "*O   R   C   A*" $line] } {
               set NameQM "ORCA"
               break
        }
    }
    close $inFile

    # Print message if wrong QM software was selected
    if { $NameQM ne $qmSoft && ( $NameQM eq "Gaussian" || $NameQM eq "ORCA" ) } {
       tk_messageBox -type ok -icon warning -message "Action halted on error!" -detail "Please change the QM selector in the main tab \
       according to the QM software you used. It seems the following file was generated by $NameQM: ${outFile}"
       set exitVar 1
    # Print message if file was not recognised
       } elseif { ! ($NameQM eq "Gaussian" || $NameQM eq "ORCA") } {
       tk_messageBox -type ok -icon warning -message "Action halted on error!" -detail "The following QM output file was not recognized: ${outFile}"
       set exitVar 1
    }
   # Return a variable that can be used to exit the current procedure
   return $exitVar

}
#======================================================
proc ::ForceFieldToolKit::SharedFcns::resetInputDefaultsESP {} {
    # resets the Input File Settings of ESP charge optiization to default values
    set ::ForceFieldToolKit::ChargeOpt::ESP::ihfree 1
    set ::ForceFieldToolKit::ChargeOpt::ESP::qwt 0.0005
    set ::ForceFieldToolKit::ChargeOpt::ESP::iqopt 2
}

#======================================================
proc ::ForceFieldToolKit::SharedFcns::checkElementPDB {} {
    # NEW 02/01/2019: Check element is defined

    set flag 0
    foreach auxAtom [[atomselect top all] get index] {
        set at [atomselect top "index $auxAtom"]
        if { [$at get element] == "X" && ![::ForceFieldToolKit::SharedFcns::LonePair::isLP $auxAtom]} {
             ::TopoTools::guessatomdata $at element mass
	     set flag 1
        }
    }
    if { $flag } {
        tk_messageBox -type ok -icon warning -message "Warning!" -detail "Some atom in your pdb file was missing the element attribute. \
        The missing elements were guessed. However, check the loaded structure or the generated output files for errors."
    }

}
#======================================================

