proc vmd_draw_arrow {mol start end col} { 
    # an arrow is made of a cylinder and a cone 
	draw color $col
    set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]] 
    graphics $mol cylinder $start $middle radius 0.15 
    graphics $mol cone $middle $end radius 0.25 
} 

proc angle_PAIM {{ mol }} {
	set span 4
	set F [ molinfo $mol get numframes]
	
	# old Davi's selection a/c' d/f' (plain wrong) set CM_sel_list [list " alpha and resid 517 264 " " alpha and resid 515 265"  " alpha and resid 267 514" " alpha and resid 512 269" " alpha and resid 270 510" " alpha and resid 508 272" " alpha and resid 274 507" " alpha and resid 505 276" " alpha and resid 277 503" " alpha and resid 501 279" " alpha and resid 281 500" " alpha and resid 498 283" " alpha and resid 284 496" " alpha and resid 494 286" " alpha and resid 288 493" " alpha and resid 491 290" " alpha and resid 291 489" " alpha and resid 487 293" " alpha and resid 295 486" " alpha and resid 484 297" " alpha and resid 298 482" " alpha and resid 480 300" " alpha and resid 302 479" " alpha and resid 477 304" " alpha and resid 305 475" " alpha and resid 473 307" " alpha and resid 309 472" " alpha and resid 470 311" " alpha and resid 312 468" " alpha and resid 466 314" " alpha and resid 316 465" " alpha and resid 463 318" " alpha and resid 319 461" " alpha and resid 459 321" " alpha and resid 323 458" " alpha and resid 456 325" " alpha and resid 326 454" " alpha and resid 452 328" " alpha and resid 330 451" " alpha and resid 449 332" " alpha and resid 333 447" " alpha and resid 445 335" " alpha and resid 337 444" " alpha and resid 442 339" " alpha and resid 340 440" " alpha and resid 438 342" " alpha and resid 344 437" " alpha and resid 435 346" " alpha and resid 347 433" " alpha and resid 431 349" " alpha and resid 351 430" " alpha and resid 428 353" " alpha and resid 354 426" " alpha and resid 424 356" " alpha and resid 358 423" " alpha and resid 421 360" " alpha and resid 361 419" " alpha and resid 417 363" " alpha and resid 365 416" " alpha and resid 414 367" " alpha and resid 368 412" " alpha and resid 410 370" " alpha and resid 372 409" " alpha and resid 407 374" " alpha and resid 375 405" " alpha and resid 403 377" " alpha and resid 379 402" " alpha and resid 400 381" " alpha and resid 382 398" " alpha and resid 396 384" " alpha and resid 386 395" " alpha and resid 393 388" " alpha and resid 389 391" ]
	#set CM_sel_list [ list " alpha and resid 265 517" " alpha and resid 267 514" " alpha and resid 269 512" " alpha and resid 272 510" " alpha and resid 274 507" " alpha and resid 276 505" " alpha and resid 279 503" " alpha and resid 281 500" " alpha and resid 283 498" " alpha and resid 286 496" " alpha and resid 288 493" " alpha and resid 290 491" " alpha and resid 293 489" " alpha and resid 295 486" " alpha and resid 297 484" " alpha and resid 300 482" " alpha and resid 302 479" " alpha and resid 304 477" " alpha and resid 307 475" " alpha and resid 309 472" " alpha and resid 311 470" " alpha and resid 314 468" " alpha and resid 316 465" " alpha and resid 318 463" " alpha and resid 321 461" " alpha and resid 323 458" " alpha and resid 325 456" " alpha and resid 328 454" " alpha and resid 330 451" " alpha and resid 332 449" " alpha and resid 335 447" " alpha and resid 337 444" " alpha and resid 339 442" " alpha and resid 342 440" " alpha and resid 344 437" " alpha and resid 346 435" " alpha and resid 349 433" " alpha and resid 351 430" " alpha and resid 353 428" " alpha and resid 356 426" " alpha and resid 358 423" " alpha and resid 360 421" " alpha and resid 363 419" " alpha and resid 365 416" " alpha and resid 367 414" " alpha and resid 370 412" " alpha and resid 372 409" " alpha and resid 374 407" " alpha and resid 377 405" " alpha and resid 379 402" " alpha and resid 381 400" " alpha and resid 384 398" " alpha and resid 386 395" " alpha and resid 388 393" ]
	#set CM_sel_list [ list " alpha and resid 265 517" " alpha and resid 272 510" " alpha and resid 279 503" " alpha and resid 286 496" " alpha and resid 293 489" " alpha and resid 300 482" " alpha and resid 307 475" " alpha and resid 314 468" " alpha and resid 321 461" " alpha and resid 328 454" " alpha and resid 335 447" " alpha and resid 342 440" " alpha and resid 349 433" " alpha and resid 356 426" " alpha and resid 363 419" " alpha and resid 370 412" " alpha and resid 377 405" " alpha and resid 384 398"]
	#set CM_sel_list [ list " alpha and resid 263 519" " alpha and resid 264 518" " alpha and resid 265 517" " alpha and resid 266 516" " alpha and resid 267 515" " alpha and resid 268 514" " alpha and resid 269 513" " alpha and resid 270 512" " alpha and resid 271 511" " alpha and resid 272 510" " alpha and resid 273 509" " alpha and resid 274 508" " alpha and resid 275 507" " alpha and resid 276 506" " alpha and resid 277 505" " alpha and resid 278 504" " alpha and resid 279 503" " alpha and resid 280 502" " alpha and resid 281 501" " alpha and resid 282 500" " alpha and resid 283 499" " alpha and resid 284 498" " alpha and resid 285 497" " alpha and resid 286 496" " alpha and resid 287 495" " alpha and resid 288 494" " alpha and resid 289 493" " alpha and resid 290 492" " alpha and resid 291 491" " alpha and resid 292 490" " alpha and resid 293 489" " alpha and resid 294 488" " alpha and resid 295 487" " alpha and resid 296 486" " alpha and resid 297 485" " alpha and resid 298 484" " alpha and resid 299 483" " alpha and resid 300 482" " alpha and resid 301 481" " alpha and resid 302 480" " alpha and resid 303 479" " alpha and resid 304 478" " alpha and resid 305 477" " alpha and resid 306 476" " alpha and resid 307 475" " alpha and resid 308 474" " alpha and resid 309 473" " alpha and resid 310 472" " alpha and resid 311 471" " alpha and resid 312 470" " alpha and resid 313 469" " alpha and resid 314 468" " alpha and resid 315 467" " alpha and resid 316 466" " alpha and resid 317 465" " alpha and resid 318 464" " alpha and resid 319 463" " alpha and resid 320 462" " alpha and resid 321 461" " alpha and resid 322 460" " alpha and resid 323 459" " alpha and resid 324 458" " alpha and resid 325 457" " alpha and resid 326 456" " alpha and resid 327 455" " alpha and resid 328 454" " alpha and resid 329 453" " alpha and resid 330 452" " alpha and resid 331 451" " alpha and resid 332 450" " alpha and resid 333 449" " alpha and resid 334 448" " alpha and resid 335 447" " alpha and resid 336 446" " alpha and resid 337 445" " alpha and resid 338 444" " alpha and resid 339 443" " alpha and resid 340 442" " alpha and resid 341 441" " alpha and resid 342 440" " alpha and resid 343 439" " alpha and resid 344 438" " alpha and resid 345 437" " alpha and resid 346 436" " alpha and resid 347 435" " alpha and resid 348 434" " alpha and resid 349 433" " alpha and resid 350 432" " alpha and resid 351 431" " alpha and resid 352 430" " alpha and resid 353 429" " alpha and resid 354 428" " alpha and resid 355 427" " alpha and resid 356 426" " alpha and resid 357 425" " alpha and resid 358 424" " alpha and resid 359 423" " alpha and resid 360 422" " alpha and resid 361 421" " alpha and resid 362 420" " alpha and resid 363 419" " alpha and resid 364 418" " alpha and resid 365 417" " alpha and resid 366 416" " alpha and resid 367 415" " alpha and resid 368 414" " alpha and resid 369 413" " alpha and resid 370 412" " alpha and resid 371 411" " alpha and resid 372 410" " alpha and resid 373 409" " alpha and resid 374 408" " alpha and resid 375 407" " alpha and resid 376 406" " alpha and resid 377 405" " alpha and resid 378 404" " alpha and resid 379 403" " alpha and resid 380 402" " alpha and resid 381 401" " alpha and resid 382 400" " alpha and resid 383 399" " alpha and resid 384 398" " alpha and resid 385 397" " alpha and resid 386 396" " alpha and resid 387 395" " alpha and resid 388 394" " alpha and resid 389 393" " alpha and resid 390 392" ]
	#set CM_sel_num [ llength $CM_sel_list ]
	#set labels [list "515 265"  "267 514" "512 269" "270 510" "508 272" "274 507" "505 276" "277 503" "501 279" "281 500" "498 283" "284 496" "494 286" "288 493" "491 290" "291 489" "487 293" "295 486" "484 297" "298 482" "480 300" "302 479" "477 304" "305 475" "473 307" "309 472" "470 311" "312 468" "466 314" "316 465" "463 318" "319 461" "459 321" "323 458" "456 325" "326 454" "452 328" "330 451" "449 332" "333 447" "445 335" "337 444" "442 339" "340 440" "438 342" "344 437" "435 346" "347 433" "431 349" "351 430" "428 353" "354 426" "424 356" "358 423" "421 360" "361 419" "417 363" "365 416" "414 367" "368 412" "410 370" "372 409" "407 374" "375 405" "403 377" "379 402" "400 381" "382 398" "396 384" "386 395" "393 388" "389 391" ]
	#set labels [list "267-514" "269-512" "272-510" "274-507" "276-505" "279-503" "281-500" "283-498" "286-496" "288-493" "290-491" "293-489" "295-486" "297-484" "300-482" "302-479" "304-477" "307-475" "309-472" "311-470" "314-468" "316-465" "318-463" "321-461" "323-458" "325-456" "328-454" "330-451" "332-449" "335-447" "337-444" "339-442" "342-440" "344-437" "346-435" "349-433" "351-430" "353-428" "356-426" "358-423" "360-421" "363-419" "365-416" "367-414" "370-412" "372-409" "374-407" "377-405" "379-402" "381-400" "384-398" "386-395" ]
	#set labels [ list "272-510" "279-503" "286-496" "293-489" "300-482" "307-475" "314-468" "321-461" "328-454" "335-447" "342-440" "349-433" "356-426" "363-419" "370-412" "377-405" ]	
	
	set labels [ list "267-515" "268-514" "269-513" "270-512" "271-511" "272-510" "273-509" "274-508" "275-507" "276-506" "277-505" "278-504" "279-503" "280-502" "281-501" "282-500" "283-499" "284-498" "285-497" "286-496" "287-495" "288-494" "289-493" "290-492" "291-491" "292-490" "293-489" "294-488" "295-487" "296-486" "297-485" "298-484" "299-483" "300-482" "301-481" "302-480" "303-479" "304-478" "305-477" "306-476" "307-475" "308-474" "309-473" "310-472" "311-471" "312-470" "313-469" "314-468" "315-467" "316-466" "317-465" "318-464" "319-463" "320-462" "321-461" "322-460" "323-459" "324-458" "325-457" "326-456" "327-455" "328-454" "329-453" "330-452" "331-451" "332-450" "333-449" "334-448" "335-447" "336-446" "337-445" "338-444" "339-443" "340-442" "341-441" "342-440" "343-439" "344-438" "345-437" "346-436" "347-435" "348-434" "349-433" "350-432" "351-431" "352-430" "353-429" "354-428" "355-427" "356-426" "357-425" "358-424" "359-423" "360-422" "361-421" "362-420" "363-419" "364-418" "365-417" "366-416" "367-415" "368-414" "369-413" "370-412" "371-411" "372-410" "373-409" "374-408" "375-407" "376-406" "377-405" "378-404" "379-403" "380-402" "381-401" "382-400" "383-399" "384-398" "385-397" "386-396"]
	set lab_num [ llength $labels ]
	set fileout [ open mcp_angle_PAIM.txt w ]
	set line ""
	
	foreach r $labels {
        set line "$line\t$r"
    }

	puts $fileout $line
	
	for { set f 0 } { $f < $F } { incr f } {
		puts "Working on frame $f of $F with $span residues"
		set line "$f"
		
		for { set i 0 } { $i < $lab_num } { incr i } {
			##puts "\tWorking on layer [ lindex $labels  $i ]"
			set res1 [lindex [ split [ lindex $labels $i ] {-} ] 0 ]
			set res2 [lindex [ split [ lindex $labels $i ] {-} ] 1 ]
			set L1 [atomselect $mol " alpha and resid [ expr $res1 - $span ] to [ expr $res1 - 1 ] [ expr $res2 + 1] to [expr $res2 + $span ] " frame $f]
			set L2 [atomselect $mol " alpha and resid [ expr $res1 + 1 ] to [ expr $res1 + $span ] [ expr $res2 - $span] to [expr $res2 - 1 ] " frame $f]
			#set L1 [ atomselect $mol "[ lindex $CM_sel_list $i ] " frame $f]
			#set L2 [ atomselect $mol "[ lindex $CM_sel_list [ expr $i + 2 ] ] " frame $f]
			set X [lindex [ lindex [measure inertia $L1] 1 ] 0 ]
			set Xu [lindex [ lindex [measure inertia $L2] 1 ] 0 ]
			if { [vecdot $X $Xu ] > 1 } { set A 0 } else {
				set A [expr 180 * acos( [vecdot $X $Xu ]) / 3.14159265358979323846 ]
			}
			if { $A > 150 } { set A [expr 180 - $A ] }
			set line [ format "$line\t%.2f" $A ]
			$L1 delete
			$L2 delete
		}
		#puts $line
		puts $fileout "$line"
	}
}

proc draw_the_thingy { {mol top } res1 res2 span} {
	#set L1 [ atomselect $mol "alpha and resid 333 to 339 442 to 449" ]
	#set L1 $sel1
	#set L2 [ atomselect $mol "alpha and resid 339 442" ]
	#set L2 $sel2
	set sel0 [atomselect $mol "alpha and resid $res1 $res2"]
	mol color ColorID 4
	mol representation VDW 0.8 20
	mol selection "alpha and resid $res1 $res2"
	mol material Glossy
	mol addrep $mol
	
	set L1 [atomselect $mol " alpha and resid [ expr $res1 - $span ] to [ expr $res1 - 1 ] [ expr $res2 + 1] to [expr $res2 + $span ] "]
	mol color ColorID 0
	mol representation CPK 1.0 0.3 10.0 10.0 
	mol selection " alpha and resid [ expr $res1 - $span ] to [ expr $res1 - 1 ] [ expr $res2 + 1] to [expr $res2 + $span ] "
	mol material Glossy
	mol addrep $mol
	set L2 [atomselect $mol " alpha and resid [ expr $res1 + 1 ] to [ expr $res1 + $span ] [ expr $res2 - $span] to [expr $res2 - 1 ] "]
	mol color ColorID 1
	mol representation CPK 1.0 0.3 10.0 10.0 
	mol selection " alpha and resid [ expr $res1 + 1 ] to [ expr $res1 + $span ] [ expr $res2 - $span] to [expr $res2 - 1 ] "
	mol material Glossy
	mol addrep $mol
	set S0 [ lindex [measure inertia $sel0] 0]
	set S [ lindex [measure inertia $L1] 0 ]
	set Su [ lindex [measure inertia $L2] 0 ]
	set X [lindex [ lindex [measure inertia $L1] 1 ] 0 ]
	set Xu [lindex [ lindex [measure inertia $L2] 1 ] 0 ]
	vmd_draw_arrow top $S0 [ vecadd $S0 [ vecscale 10 $X ]] blue
	vmd_draw_arrow top $S0 [ vecadd $S0 [ vecscale 10 $Xu ]] red
	#vmd_draw_arrow top $S [ vecadd $S [ vecscale 10 $X ]] blue
	#vmd_draw_arrow top $Su [ vecadd $Su [ vecscale 10 $Xu ]] red
	if { [vecdot $X $Xu ] > 1 } { set A 0 } else {
				set A [expr 180 * acos( [vecdot $X $Xu ]) / 3.14159265358979323846 ]
			}
			if { $A > 150 } { set A [expr 180 - $A ] }
		
	puts $A
	}
