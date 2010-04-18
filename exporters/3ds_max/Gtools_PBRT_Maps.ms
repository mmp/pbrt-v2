-- // =======================================================================================================================================================
-- // Constant

plugin TextureMap PBRT_ConstantMap name:"PBRT Constant Map" classID:#(0x69313320, 0x34cb02ca) extends:Bitmaptexture replaceUI:true version:1
( 
	parameters main rollout:params
	(
		color type:#color default:[100,100,100] ui:m_Color
	)

	rollout params "Constant Parameters"
	(
		colorpicker m_Color "Color: " align:#left
	)
)

-- // =======================================================================================================================================================
-- // Checkerboard

plugin TextureMap PBRT_CheckboardMap name:"PBRT Checkboard Map" classID:#(0x21c0d105, 0x6ee5c2e7) extends:Bitmaptexture replaceUI:true version:1
( 
	parameters main rollout:params
	(
		aamode type:#integer default:1 ui:m_aamode
		color1 type:#color default:[0,0,0] ui:m_Color1
		color2 type:#color default:[255,255,255] ui:m_Color2
	)

	parameters main rollout:paramsCoordinates
	(
		uv type:#integer default:1 ui:m_uv
		uscale type:#float default:1 ui:m_uscale
		vscale type:#float default:1 ui:m_vscale
		udelta type:#float default:0 ui:m_udelta
		vdelta type:#float default:0 ui:m_vdelta
	)
	
	rollout params "Checkboard Parameters"
	(
		dropdownlist m_aamode "AA Mode:" items:#("Closedform", "SuperSample", "None") align:#l width:200
		colorpicker m_Color1 "Color: " align:#left
		colorpicker m_Color2 "Color: " align:#left
	)

	rollout paramsCoordinates "Coordinates"
	(
		dropdownlist m_uv "Mapping:" items:#("UV", "Spherical", "Cylindrical", "Planar") align:#right width:200

		spinner m_uscale "U Scale:" range:[-100, 100, 1] width:100 align:#right
		spinner m_vscale "V Scale:" range:[-100, 100, 1] width:100 align:#right
		
		spinner m_udelta "U Delta:" range:[-100, 100, 0] width:100 align:#right
		spinner m_vdelta "V Delta:" range:[-100, 100, 0] width:100 align:#right
	)
	
	on create do
	(
	)
)

-- // =======================================================================================================================================================
-- // ImageMap

plugin TextureMap PBRT_BitMap name:"PBRT Bitmap" classID:#(0x6d0f0ce6, 0x28c251dc) extends:Bitmaptexture replaceUI:true version:1
( 
	local loaded = false
	
	parameters main rollout:params
	(
		
		pbrt_bitmap type:#bitmap ui:Imap_bitmap
		wrap type:#integer default:1 ui:Imap_wrap
		tri type:#boolean default:false ui:Imap_Tri_chk
		anist type:#float deafult:8.0 	ui:Imap_Anis_spn
		gamma type:#float deafult:1.0 ui:Imap_Gamma_spn 

--		_scale			type:#float			deafult:1.0		ui:Imap_Scale_spn
		
		on pbrt_bitmap set val do (if loaded then delegate.bitmap = val)


	)

	parameters main rollout:paramsCoordinates
	(
		uv type:#integer default:1 ui:m_uv
		uscale type:#float default:1 ui:m_uscale
		vscale type:#float default:1 ui:m_vscale
		udelta type:#float default:0 ui:m_udelta
		vdelta type:#float default:0 ui:m_vdelta
		
		on uscale set val do delegate.coords.u_tiling = val
		on vscale set val do delegate.coords.v_tiling = val	
		on udelta set val do delegate.coords.u_offset = val
		on vdelta set val do delegate.coords.v_offset = val				
	)

-- ImageTexture(TextureMapping2D *m, const string &filename, bool doTrilinear, float maxAniso, ImageWrap wrapMode, float scale, float gamma)		
	rollout params "ImageMap Parameters"
	(
		label 			lb_filename "Bitmap:" pos:[8,10]
		button 			Imap_bitmap "None" pos:[50,8] width:265
		label 			Imap_Wrap_lbl "Wrap:" pos:[8,40]
		dropdownlist Imap_wrap "" items:#("Repeat", "Black", "Clamp") pos:[40,38] width:70
		checkbox		Imap_Tri_chk "Use Trilinear Interpolation" checked:false pos:[140,42] width:200 
		spinner			Imap_Anis_spn "Max Anisotropy:" range:[1.0, 512.0, 8.0] pos:[60,70] width:80 
		spinner			Imap_Gamma_spn "Gamma:" pos:[30,95] width:80	
		button			Imap_seeBitmap_btn "View Bitmap"  pos:[200,90] width:70 
--		spinner			Imap_Scale_spn "Scale:" range:[0.0, 10.0, 1.0] pos:[32,117] width:80

		
		on params open do (
/*	*/			
					try(
						Imap_bitmap.text = pbrt_bitmap.filename
						Imap_bitmap.tooltip = pbrt_bitmap.filename
						)catch()
		
		)
	
		on Imap_seeBitmap_btn pressed do (
			try(display PBRT_BitMap)catch()
		)
		
		on Imap_bitmap pressed do (
/**/
					loaded = true
					try(
					pbrt_bitmap = selectBitMap()
					Imap_bitmap.text = pbrt_bitmap.filename
					Imap_bitmap.tooltip = pbrt_bitmap.filename
						)catch()	

		)
	
	)
	
	on load do (
--		print "load"
		delegate.bitmap	
	)
	
	on create do (
		del = delegate
		-- for some reason the values dont get their default value so set it up here.
		this.anist = 8.0
		this.gamma = 1.0
	)
	
	on clone old do (
		del = delegate
--		print "clone"
	)
	
	
	rollout paramsCoordinates "Coordinates"
	(
		dropdownlist m_uv "Mapping:" items:#("UV", "Spherical", "Cylindrical", "Planar") align:#right width:200

		spinner m_uscale "U Scale:" range:[-100, 100, 1] width:100 align:#right
		spinner m_vscale "V Scale:" range:[-100, 100, 1] width:100 align:#right
		
		spinner m_udelta "U Delta:" range:[-100, 100, 0] width:100 align:#right
		spinner m_vdelta "V Delta:" range:[-100, 100, 0] width:100 align:#right
		
		on params open do (
			
					try(
						m_uscale.value = pbrt_bitmap.coords.u_tiling
						m_vscale.value = pbrt_bitmap.coords.v_tiling
						m_udelta.value = pbrt_bitmap.coords.u_offset
						m_vdelta.value = pbrt_bitmap.coords.v_offset
						)catch()
		
		)
	)
)