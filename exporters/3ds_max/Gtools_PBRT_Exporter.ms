/*
Started to work on it. Nov 17 2009

PBRT 2 Exporter by Guillermo M. Leal LLaguno.  
g.llaguno@evvisual.com
www.evvisual.com

Script to export Max geometry, lights etc to pbrt v2 render, http://www.pbrt.org/  , http://github.com/mmp/pbrt-v2

You can use the script anyway you like, but please mention my name/website if you use it. And if you do something cool I would love to see. And last if you made 
	Changes/modification please email it to me.
	
I created this script so I can create the pbrt v2 book cover, having said that it was done in a rush so it’s incomplete and not very clean but it works.

Reference:
	- pbrt book.
	- Matt Pharr .  Many question i had to ask... thanks.
	- Blender luxrender python exporter.
	- Luxrender wiki and forum.

The script is slow.. Specially on large files, but i had no time to optimize, I think part is that i add a modifier to turn to mesh and 
Then remove it to every object. I was having a lot of problems with instances and for some reason doing this solve the issue.

The render button is disable by default, the reason is that if you render from within max, max will be frozen until the render finishes plus you won’t see all
	Information pbrt give you.  It’s better to export it and then open the command line to run it, like  c:\pbrt test.pbrt


Current futures/implementations

-	Mesh objects with uv's, if the object is hidden it’s not exported.
-	Lights:  
			- Point Lights represented by omni lights, it uses the position, color and multiplier of the omni light
			- Spot Lights represented by Target Spot Lights, it uses the position, color, hotspot, falloff and multiplier of the Spot Target Lights
			- Distant Lights represented by target Directional lights, it uses the position, color and multiplier of the omni light
-  Area Lights from geometry it uses the PBRT_Light material, for rectangles its fast for complex geometry it’s too slow.
- Environment Color
- Environment Map

- Materials:
	- pbrt only allows exr and tga files, to use a bitmap you must use PBRT Bitmap map.
	- If you want opacity maps, the need to be a black/white or (gray tones) tga or exr , that is not an alpha channel.
	- Almost all materials support opacity map and bumpmap.
	- the amounts for the parameters are real low numbers, not like max let say in max you my use 30 for bump, in pbrt it might be .001

- PBRT Materials currently supported:
			- Matte
			- Plastic
			- Light Material (area)
			- Constant
			- Uber
			- Glass
			- Translucent
			- Metal
			- Substrate
			
- Textures: 	
			- PBRT Bitmap
			- PBRT CheckerBoard
			- PBRT Constant
			
- Cameras:
			- Perspective
- Instances


The render settings are saved in a point helper, named PBRT_Helper if it’s deleted it will be recreated next time you run the script.
	This part of the script its not finish there are 2 functions PBRT_SaveProperties and PBRT_GetProperties they need to updated so all settings can be saved
	At the moment only the render tab are saved.
	
Brief description of the script parameters.

System Tab:
	- Save Path: where you want the exported .pbrt file to be and name, don’t include extension.
	- PBRT exe: location of pbrt.exe, this can be hardcode it in the variable PBRT_EXE_PATH
	- Threads: number of threads used by PBRT
	- use Instances: in case you don’t want to use?
	- disable maps. The materials will get exported but no maps, only the diffuse color, this is useful sometimes as the textures are districting to detect some issue.
	- Accelerator type.
	
Output Tab:
	- width:
	- height:
	- crop: use max crop rectangle, (handy to have)
	
	
	-- ISUES
	- I need to add check for undefined materials

*/





(		
--  ***************************************************************
Global PBRTMainFloater
try(closeRolloutFloater PBRTMainFloater)catch()
		
--Local Viewer = "\"" + "C:\\Users\\Memo\\Desktop\\MisChingaderas\\Gtools_C++\\Libs\\OpenEXR\\Deploy\\bin_VC2008\\Release\\"
Local Viewer = "\"" + "C:\\Program Files (x86)\\exr\\"
	
-- GUI
Global PBRT_System
Global PBRT_Output
Global PBRT_Render
Global PBRT_Camera
Global PBRT_Tools
--
Local fname = "pbrt_test"		
Local SAVEPATH =  (GetDir #export) + "\\" + fname
-- Local PBRT_EXE_PATH = "C:\\Users\\Memo\\Desktop\\MisChingaderas\\Gtools_C++\\3D\\Raytracing\\pbrt\\PBRT_2_Bin\\"
--Local PBRT_EXE_PATH = "O:\\"
	
Local ResRatio
Local CurrentCamera = undefined
Local doneMats = #() -- array to hold the materials that have been already processed
Local TextArray = #()	-- array to hold the proccesed textures
Local GlobalMaterial = undefined -- material to use for overraide
--	Local out_file
ver = "1.00"

	
fn CopyMaps thepath = (
	
		theMaps = #()
	
		TheObjects = for o in objects collect o
		
		-- look for all maps
		for o = 1 to TheObjects.count do (
			
			if TheObjects[o].material != undefined then (-- check if it has a material
					Local mat = TheObjects[o].material
					-- matte material
					if classof mat == PBRT_Matte then (
						
							-- KdTexture
							-- check if it has a texture and its a pbrt texture
							if mat.KdTexture != undefined and classof mat.KdTexture == PBRT_BitMap then (
								local theBitmap = mat.KdTexture.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)
						
							-- alphaImage
							-- check if it has a texture and its a pbrt texture
							if mat.alphaImage != undefined and classof mat.alphaImage == PBRT_BitMap then (
								local theBitmap = mat.alphaImage.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)
						
							-- sigmaMap
							-- check if it has a texture and its a pbrt texture
							if mat.sigmaMap != undefined and classof mat.sigmaMap == PBRT_BitMap then (
								local theBitmap = mat.sigmaMap.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)
						
							-- bumpMap
							-- check if it has a texture and its a pbrt texture
							if mat.bumpMap != undefined and classof mat.bumpMap == PBRT_BitMap then (
								local theBitmap = mat.bumpMap.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)								
						) -- end matte material
						
					-- plastic material
					if classof mat == PBRT_Plastic then (
						
							-- KdTexture
							-- check if it has a texture and its a pbrt texture
							if mat.KdTexture != undefined and classof mat.KdTexture == PBRT_BitMap then (
								local theBitmap = mat.KdTexture.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)
						
							-- alphaImage
							-- check if it has a texture and its a pbrt texture
							if mat.alphaImage != undefined and classof mat.alphaImage == PBRT_BitMap then (
								local theBitmap = mat.alphaImage.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)
						
							-- KsTexture
							-- check if it has a texture and its a pbrt texture
							if mat.KsTexture != undefined and classof mat.KsTexture == PBRT_BitMap then (
								local theBitmap = mat.KsTexture.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)
						
							-- RoughnessTexture
							-- check if it has a texture and its a pbrt texture
							if mat.RoughnessTexture != undefined and classof mat.RoughnessTexture == PBRT_BitMap then (
								local theBitmap = mat.RoughnessTexture.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)
								
							-- bumpMap
							-- check if it has a texture and its a pbrt texture
							if mat.bumpMap != undefined and classof mat.bumpMap == PBRT_BitMap then (
								local theBitmap = mat.bumpMap.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)								
						) -- end plastic material
						
					-- Uber material
					if classof mat == PBRT_Uber then (
						
							-- KdTexture
							-- check if it has a texture and its a pbrt texture
							if mat.KdTexture != undefined and classof mat.KdTexture == PBRT_BitMap then (
								local theBitmap = mat.KdTexture.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)
						
							-- KsTexture
							-- check if it has a texture and its a pbrt texture
							if mat.KsTexture != undefined and classof mat.KsTexture == PBRT_BitMap then (
								local theBitmap = mat.KsTexture.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)	
								
							-- opacityTexture
							-- check if it has a texture and its a pbrt texture
							if mat.opacityTexture != undefined and classof mat.opacityTexture == PBRT_BitMap then (
								local theBitmap = mat.opacityTexture.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)
						
							-- RoughnessTexture
							-- check if it has a texture and its a pbrt texture
							if mat.RoughnessTexture != undefined and classof mat.RoughnessTexture == PBRT_BitMap then (
								local theBitmap = mat.RoughnessTexture.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)
								
							-- bumpMap
							-- check if it has a texture and its a pbrt texture
							if mat.bumpMap != undefined and classof mat.bumpMap == PBRT_BitMap then (
								local theBitmap = mat.bumpMap.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)								
						) -- end Uber material
						
					-- Uber material
					if classof mat == PBRT_Metal then (
						

								
							-- eta
							if mat.eta != undefined then (
								local thefile = mat.eta
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if thefile != undefined and (finditem theMaps  thefile) == 0 then append theMaps thefile
								)

							-- k
							if mat.k != undefined then (
								local thefile = mat.k
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if thefile != undefined and (finditem theMaps  thefile) == 0 then append theMaps thefile
								)
								
							-- RoughnessTexture
							-- check if it has a texture and its a pbrt texture
							if mat.RoughnessTexture != undefined and classof mat.RoughnessTexture == PBRT_BitMap then (
								local theBitmap = mat.RoughnessTexture.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)
								
							-- bumpMap
							-- check if it has a texture and its a pbrt texture
							if mat.bumpMap != undefined and classof mat.bumpMap == PBRT_BitMap then (
								local theBitmap = mat.bumpMap.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)								
						) -- end Uber material						
					
					-- Translucent material
					if classof mat == PBRT_Translucent then (
						
							-- KdTexture
							-- check if it has a texture and its a pbrt texture
							if mat.KdTexture != undefined and classof mat.KdTexture == PBRT_BitMap then (
								local theBitmap = mat.KdTexture.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)
						
							-- KsTexture
							-- check if it has a texture and its a pbrt texture
							if mat.KsTexture != undefined and classof mat.KsTexture == PBRT_BitMap then (
								local theBitmap = mat.KsTexture.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)	
								
							-- alphaImage
							-- check if it has a texture and its a pbrt texture
							if mat.alphaImage != undefined and classof mat.alphaImage == PBRT_BitMap then (
								local theBitmap = mat.alphaImage.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)
						
							-- RoughnessTexture
							-- check if it has a texture and its a pbrt texture
							if mat.RoughnessTexture != undefined and classof mat.RoughnessTexture == PBRT_BitMap then (
								local theBitmap = mat.RoughnessTexture.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)
								
							-- reflectTexture
							-- check if it has a texture and its a pbrt texture
							if mat.reflectTexture != undefined and classof mat.reflectTexture == PBRT_BitMap then (
								local theBitmap = mat.reflectTexture.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)
								
							-- transmitTexture
							-- check if it has a texture and its a pbrt texture
							if mat.transmitTexture != undefined and classof mat.transmitTexture == PBRT_BitMap then (
								local theBitmap = mat.transmitTexture.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)
								
							-- bumpMap
							-- check if it has a texture and its a pbrt texture
							if mat.bumpMap != undefined and classof mat.bumpMap == PBRT_BitMap then (
								local theBitmap = mat.bumpMap.pbrt_bitmap.filename
								-- if the bitmap its not undefined and its not in the maps array add it to the array
								if theBitmap != undefined and (finditem theMaps  theBitmap) == 0 then append theMaps theBitmap
								)								
						) -- end Translucent material
						
				) -- end if it has material
		) -- end object loop

/*		*/
		-- copy the maps
		for c = 1 to theMaps.count do (
			theFile = thepath + (filenameFromPath theMaps[c])
			copyfile theMaps[c] theFile -- copy the files
			)

--		print theMaps
	)
	
------------------------------------------------
-- this is not my function but i cant remember where i got i from,  sorry. thanks who every did it.
fn clampRange 
	stVal:0.0  			--Start value of in comming range.
	edVal:100.0 		--End value of incomming range.
	nStVal:0.0 			--Start value of out going range.
	nEdVal:1.0 			--End value of out going range.
	val:0.0 			--In value that will be converted.
	clamp:false			--If true will return a value that is clamped in the out going range. 
	=
	(
		local outVal
		
		--Get the percent of the inVal
		inRange=(edVal as float)-(stVal as float)
		inVal=(val as float)-(stVal as float)
		valPercent=inVal/inRange
		
		--Get the value of the based on the new range.
		nRange=(nEdVal as float)-(nStVal as float)
		
		--Invert the math if the nStVal> nEdVal so that values ranges can be inverted.
		if nStVal<=nEdVal then
		(
			outVal=(nRange*valPercent)-(abs nStVal)
		)else
		(
			outVal=(nRange*valPercent)+(abs nStVal)
		)
		
		--Clamp out going value to out range if clamp is true.
		if clamp then
		(
			if outVal<nStVal then (return nStVal)
			if outVal>nEdval then (return nEdVal)
		)
		return outVal
	)

fn PbrtColor c =
	(
		return ( ((c.r / 255)) as string + " " + ((c.g / 255)) as string + " " + ((c.b / 255)) as string )
	)

fn PbrtAddBack str =
	(
		if (str == undefined) do 
			str = ""
		
	--	str = substituteString str "\\" "\\\\"
		str = substituteString str "\\" "/"
		
		return ("\"" + str + "\"") 
	)

-- Camera
fn Pbrt_ExpCamera cam &out_file = (
		
		-- we need to change it to  trget camera to calculate the distance
--		cam.type = #free

	
		cam_pos 	= cam.transform.translationpart
--		cam_pos 	= in coordsys local cam.pos
		cam_dir		= cam.dir
		cam_look 	= cam_pos - cam_dir
		cam_up		= cam.transform.row2
	
		camRadius = PBRT_Camera.PBRT_Camera_LensR_spn.value
		camDistance = PBRT_Camera.PBRT_Camera_FDist_spn.value
	
		format "# Camera %\n" 	cam.name to:out_file
		format "LookAt % % %\n" 	cam_pos[1] cam_pos[2] cam_pos[3] to:out_file
		format "\t\t % % %\n" 	cam_look[1] cam_look[2] cam_look[3] to:out_file
		format "\t\t % % %\n\n" 	cam_up[1] cam_up[2] cam_up[3] to:out_file
/*	
		camera_transform = inverse (GetViewTM())
		position = camera_transform.position
		up_vector = normalize (camera_transform.row2)
		look_vector = camera_transform.position + (normalize (-camera_transform.row3))			
			
		format "LookAt % % %\n"position.x position.y position.z  to:out_file		
		format "% % % \n" look_vector.x look_vector.y look_vector.z  to:out_file
		format "% % % \n"up_vector.x up_vector.y up_vector.z to:out_file		
			

			format "Transform [\n" to:out_file 
			format "% % % 0 \n" cam.transform[1][1] cam.transform[1][2] cam.transform[1][3]  to:out_file 
			format "% % % 0 \n" cam.transform[2][1] cam.transform[2][2] cam.transform[2][3]   to:out_file 
			format "% % % 0 \n" cam.transform[3][1] cam.transform[3][2] cam.transform[3][3]   to:out_file 
			format "% % % 1 \n" cam.transform[4][1] cam.transform[4][2] cam.transform[4][3]   to:out_file 
			format "]\n" to:out_file 
	
		format "\n" 	to:out_file
		*/	
		format "Camera \"perspective\" \n"  to:out_file
		format "\t \"float fov\" [%] \n"  cam.fov to:out_file
--		format "\t \"float hiter\" [%] \n"  cam.nearclip to:out_file
--		format "\t \"float yon\" [%] \n"  cam.farclip to:out_file
		format "\t \"float shutteropen\" [%] \n" PBRT_Camera.PBRT_Camera_ShOpen_spn.value  to:out_file
		format "\t \"float shutterclose\" [%] \n"  PBRT_Camera.PBRT_Camera_ShClose_spn.value to:out_file		
		format "\t \"float lensradius\" [%] \n"   camRadius to:out_file		
		format "\t \"float focaldistance\" [%] \n" camDistance  to:out_file			

		format "\n" 	to:out_file
		-- Render Settings
		format "Film \"image\" \n"   to:out_file
		format "\t \"integer xresolution\"  [%]\n"  PBRT_Output.RenderWidth_spn.value to:out_file
		format "\t \"integer yresolution\"  [%]\n"  PBRT_Output.RenderHeight_spn.value to:out_file
		
--		format "\t \"integer writefrequency\" [%] \n" PBRT_System.PBRT_WriteFrequency.value  to:out_file
		
		maxCrop = viewport.getRegionRect 1
	
		min_X = clampRange stVal:0 edVal:renderwidth nStVal:0.0 nEdVal:1.0 val:maxCrop.left 
		max_X = clampRange stVal:0 edVal:renderwidth nStVal:0.0 nEdVal:1.0 val:maxCrop.right 
		min_Y = clampRange stVal:0 edVal:renderheight nStVal:0.0 nEdVal:1.0 val:maxCrop.top 
		max_Y = clampRange stVal:0 edVal:renderheight nStVal:0.0 nEdVal:1.0 val:maxCrop.bottom 
		
		if PBRT_Output.crop_btn.checked then format "\t \"float cropwindow\"  [%  %  %  %]\n"  min_X max_X min_Y max_Y to:out_file
	
	
		Local theFileOnly = filterstring PBRT_System.PBRT_SavePath_edt.text "\\"			
		
			-- we need to check if the user provided file extension so we dont duplicat it
	Local theName = theFileOnly[theFileOnly.count] 	-- name pat of the complete string with or without extension, previously we already spli it by \
	Local remExt =  filterstring theName "."	-- split it again now by .

	-- if it has cont more then 1 it means it has an extension
	if remExt.count > 1 then theFileOnly = remExt[1] -- assign the element 1, removing the extension
		else theFileOnly = remExt[remExt.count] 		-- the user did not provide the extension

	
		Local TheTexFile = theFileOnly + PBRT_System.PBRT_Ext_ddl.selected
		format "\t \"string filename\"  [\"%\"]\n" TheTexFile  to:out_file  --RenderFileName RenderFileName = fname + ".exr"
		
		format "# End Camera %\n" 	cam.name to:out_file
	)

fn Pbrt_Enviorment &out_file = (
	
	envColor = (backgroundColor / 255) --* PBRT_Render.envMult.value
	
	format "\n# Enviorment \n\n" to:out_file
	
	-- if there is a map use it, if not use the color
	if useEnvironmentMap == on and EnvironmentMap != undefined then (
		theMap = EnvironmentMap.bitmap.filename
		theMapSplit = filterstring theMap "\\" -- split the filename by the \
		newPath = ""

		for o = 1 to theMapSplit.count do (
			if o != theMapSplit.count then newPath = newPath + theMapSplit[o] + "\\\\" -- now add double \\ so pbrt can read it, to all words expet the lastone
				else newPath = newPath + theMapSplit[o] -- now add the actual filename
			)
		theresult = clampRange stVal:0.0 edVal:1.0 nStVal:0.0 nEdVal:360.0 val:EnvironmentMap.coords.u_offset
		theresult += 90
		format "AttributeBegin \n" to:out_file	
		format "Rotate % 0 0 1\n"  theresult to:out_file
		format "LightSource \"infinite\" \n" to:out_file	
		format "\"integer nsamples\" [%]\n" PBRT_Render.envSamp.value  to:out_file
		format "\"color scale\" [% % %]\n" PBRT_Render.envMult.value PBRT_Render.envMult.value PBRT_Render.envMult.value  to:out_file	
		format "\t\"string mapname\" [\"%\"]\n"  newPath to:out_file
	)else (
			
		format "AttributeBegin \n" to:out_file	
		format "LightSource \"infinite\" \n" to:out_file	
		format "\t\"color L\" [% % %]\n"  envColor.r envColor.g envColor.b to:out_file
		format "\"integer nsamples\" [%]\n" PBRT_Render.envSamp.value  to:out_file
		format "\"color scale\" [% % %]\n" PBRT_Render.envMult.value PBRT_Render.envMult.value PBRT_Render.envMult.value  to:out_file	
		)
--	format "\t\"float scale\" [1.000000]\n" to:out_file
	format "AttributeEnd \n\n" to:out_file
	
)

-- filter type
fn Pbrt_PixelFilter filType &out_file = (	

		format "\n" 	to:out_file
		-- get settings
		if filType == "box" then
			(
				format "PixelFilter \"box\" \n" to:out_file
				format "\t \"float xwidth\" [%] \n" PBRT_Render.Render_Filter_Box_X.value to:out_file
				format "\t \"float ywidth\" [%] \n" PBRT_Render.Render_Filter_Box_Y.value to:out_file
			)		
		else if filType == "gaussian" then
			(
				format "PixelFilter \"gaussian\" \n" to:out_file
				format "\t \"float xwidth\" [%] \n" PBRT_Render.Render_Filter_Gauss_X.value to:out_file
				format "\t \"float ywidth\" [%] \n" PBRT_Render.Render_Filter_Gauss_Y.value to:out_file
--				format "\t \"alpha\" [%] \n" PBRT_Render.Render_Filter_Gauss_Alpah.value to:out_file
			)
			
		else if filType == "mitchell" then
			(
				local B = 1.0 - PBRT_Render.Render_Filter_Mitch_Sharp.value
				local C = PBRT_Render.Render_Filter_Mitch_Sharp.value * .5
				
				format "PixelFilter \"mitchell\" \n" to:out_file
				format "\t \"float xwidth\" [%] \n" PBRT_Render.Render_Filter_Mitch_X.value to:out_file
				format "\t \"float ywidth\" [%] \n" PBRT_Render.Render_Filter_Mitch_Y.value to:out_file
				format "\t \"float B\" [%] \n" B to:out_file
				format "\t \"float C\" [%] \n" C to:out_file
			)
			
		else if filType == "sinc" then
			(
				format "PixelFilter \"sinc\" \n" to:out_file
				format "\t \"float xwidth\" [%] \n" PBRT_Render.Render_Filter_Sinc_X.value to:out_file
				format "\t \"float ywidth\" [%] \n" PBRT_Render.Render_Filter_Sinc_Y.value to:out_file
				format "\t \"tau\" [%] \n" PBRT_Render.Render_Filter_Sinc_Tau.value to:out_file
			)
			
		else if filType == "triangle" then
			(
				format "PixelFilter \"triangle\" \n" to:out_file
				format "\t \"float xwidth\" [%] \n" PBRT_Render.Render_Filter_Tri_X.value to:out_file
				format "\t \"float ywidth\" [%] \n" PBRT_Render.Render_Filter_Tri_Y.value to:out_file
			)		
	)	
	
-- Sampler type
fn Pbrt_Sampler samType &out_file = (	

	format "\n" 	to:out_file

	-- get settings
	if samType == 1 then -- "lowdiscrepancy" then
		(
			format "Sampler \"lowdiscrepancy\" \n" to:out_file
			format "\t \"integer pixelsamples\" [%] \n" PBRT_Render.Render_Sampler_spn.value to:out_file
		)		
	else if samType == 2 then --"random" then
		(
			format "Sampler \"random\" \n" to:out_file
			format "\t \"integer pixelsamples\" [%] \n" PBRT_Render.Render_Sampler_spn.value to:out_file
		)
	)
		
-- Integrator type
fn Pbrt_SurfaceIntegrator inteType &out_file = (

	format "\n" 	to:out_file
		
	-- get settings
	if inteType == 1 then -- Ambient Occlusion
		(
			format "SurfaceIntegrator \"ambientocclusion\" \n" to:out_file
			format "\t \"integer nsamples\" [%] \n" PBRT_Render.Render_Integ_AO_Samples_spn.value to:out_file
			format "\t \"float maxdist\" [%] \n" PBRT_Render.Render_Integ_AO_MaxDist_spn.value to:out_file
		)	
	else if inteType == 2 then -- Diffuse Seprt
		(
			format "SurfaceIntegrator \"diffuseprt\" \n" to:out_file
			format "\t \"integer lmax\" [%] \n" PBRT_Render.Render_Integ_Diffuse_Max_spn.value to:out_file
			format "\t \"integer nsamples\" [%] \n" PBRT_Render.Render_Integ_Diffuse_Ns_spn.value to:out_file
		)	
	else if inteType == 3 then -- Dipole Sub Surface
		(
			format "SurfaceIntegrator \"dipolesubsurface\" \n" to:out_file
			format "\t \"integer maxdepth\" [%] \n" PBRT_Render.Render_Integ_Dipole_MaxDepth_spn.value to:out_file
			format "\t \"float maxerror\" [%] \n" PBRT_Render.Render_Integ_Dipole_MaxError_spn.value to:out_file
			format "\t \"float minsampledistance\" [%] \n" PBRT_Render.Render_Integ_Dipole_MinDist_spn.value to:out_file
		)	
	else if inteType == 4 then -- Direct Lighting
		(
			format "SurfaceIntegrator \"directlighting\" \n" to:out_file
			format "\t \"integer maxdepth\" [%] \n" PBRT_Render.Render_Integ_Direct_Depth_spn.value to:out_file
		)	
	else if inteType == 5 then -- Emission
		(
			format "SurfaceIntegrator \"emission\" \n" to:out_file
			format "\t \"float stepsize\" [%] \n" PBRT_Render.Render_Integ_Emission_Step_spn.value to:out_file
		)	
	else if inteType == 6 then -- Glossy Prt
		(
--			local Kd = PbrtColor mat.Kd
			
			format "SurfaceIntegrator \"glossyprt\" \n" to:out_file
			format "\t \"integer lmax\" [%] \n" PBRT_Render.Render_Integ_Glossy_LMax_spn.value to:out_file
			format "\t \"integer nsamples\" [%] \n" PBRT_Render.Render_Integ_Glossy_Ns_spn.value to:out_file
			format "\t \"bool dotransfer\" [\"%\"] \n" PBRT_Render.Render_Integ_Glossy_DotR_btn.state to:out_file
			format "\t \"float Kd\" [%] \n" PBRT_Render.Render_Integ_Glossy_Kd_spn.value to:out_file
			format "\t \"float Ks\" [%] \n" PBRT_Render.Render_Integ_Glossy_Ks_spn.value to:out_file
			format "\t \"float roughness\" [%] \n" PBRT_Render.Render_Integ_Glossy_Rough_spn.value to:out_file
		)	
	else if inteType == 7 then -- Instant Global Illumination
		(
--			local Kd = PbrtColor mat.Kd
			
			format "SurfaceIntegrator \"igi\" \n" to:out_file
			format "\t \"integer nlights\" [%] \n" PBRT_Render.Render_Integ_IGI_LightPath_spn.value to:out_file
			format "\t \"integer nsets\" [%] \n" PBRT_Render.Render_Integ_IGI_LightSet_spn.value to:out_file
			format "\t \"float mindist\" [%] \n" PBRT_Render.Render_Integ_IGI_MinDist_spn.value to:out_file
			format "\t \"float rrthreshold\" [%] \n" PBRT_Render.Render_Integ_IGI_RRThrs_spn.value to:out_file
			format "\t \"integer maxdepth\" [%] \n" PBRT_Render.Render_Integ_IGI_MaxDepth_spn.value to:out_file			
			format "\t \"float glimit\" [%] \n" PBRT_Render.Render_Integ_IGI_GLimit_spn.value to:out_file
			format "\t \"integer gathersamples\" [%] \n" PBRT_Render.Render_Integ_IGI_Gather_spn.value to:out_file			
		)
	else if inteType == 8 then -- Irradiance Cache
		(
			format "SurfaceIntegrator \"irradiancecache\" \n" to:out_file
			format "\t \"float minweight\" [%] \n" PBRT_Render.Render_Integ_IrCach_MinWht_spn.value  to:out_file			
			format "\t \"float minpixelspacing\" [%] \n" PBRT_Render.Render_Integ_IrCach_MinSpac_spn.value  to:out_file
			format "\t \"float maxpixelspacing\" [%] \n" PBRT_Render.Render_Integ_IrCach_MaxSpac_spn.value  to:out_file
			format "\t \"float maxangledifference\" [%] \n" PBRT_Render.Render_Integ_IrCach_MaxAngle_spn.value  to:out_file			
			format "\t \"integer axspeculardepth\" [%] \n" PBRT_Render.Render_Integ_IrCach_SpecDepth_spn.value to:out_file
			format "\t \"integer maxindirectdepth\" [%] \n" PBRT_Render.Render_Integ_IrCach_indDepth_spn.value to:out_file			
			format "\t \"integer nsamples\" [%] \n" PBRT_Render.Render_Integ_IrCach_Samples_spn.value to:out_file				
		)				
	else if inteType == 9 then -- Path Tracing
		(
			format "SurfaceIntegrator \"path\" \n" to:out_file
			format "\t \"integer maxdepth\" [%] \n" PBRT_Render.Render_Integ_Path_Depth_spn.value to:out_file
		)
	else if inteType == 10 then -- Photon Map
		(
			format "SurfaceIntegrator \"photonmap\" \n" to:out_file
			format "\t \"integer causticphotons\" [%] \n" PBRT_Render.Render_Integ_Photon_Caustic_spn.value to:out_file
			format "\t \"integer indirectphotons\" [%] \n" PBRT_Render.Render_Integ_Photon_Indirect_spn.value to:out_file			
			format "\t \"integer nused\" [%] \n" PBRT_Render.Render_Integ_Photon_nUsed_spn.value to:out_file	
			format "\t \"integer maxspeculardepth\" [%] \n" PBRT_Render.Render_Integ_Photon_SpecDepth_spn.value to:out_file			
			format "\t \"integer maxphotondepth\" [%] \n" PBRT_Render.Render_Integ_Photon_MaxPhoDepth_spn.value to:out_file	
			format "\t \"bool finalgather\" [\"%\"] \n" PBRT_Render.Render_Integ_Photon_FinalGath_btn.state to:out_file	
			format "\t \"integer finalgathersamples\" [%] \n" PBRT_Render.Render_Integ_Photon_gatherSamp_spn.value to:out_file	
			format "\t \"float maxdist\" [%] \n" PBRT_Render.Render_Integ_Photon_MaxDist_spn.value  to:out_file			
			format "\t \"float gatherangle\" [%] \n" PBRT_Render.Render_Integ_Photon_GathAng_spn.value  to:out_file			
		)
	else if inteType == 11 then -- Single
		(
			format "SurfaceIntegrator \"single\" \n" to:out_file		
			format "\t \"float stepsize\" [%] \n" PBRT_Render.Render_Integ_Single_Step_spn.value  to:out_file			
		)
	else if inteType == 12 then -- Use Probes
		(
			format "SurfaceIntegrator \"useprobes\" \n" to:out_file
		)	
	else if inteType == 13 then -- whitted
		(
			format "SurfaceIntegrator \"whitted\" \n" to:out_file		
			format "\t \"integer maxdepth\" [%] \n" PBRT_Render.Render_Integ_Whitted_MaDepth_spn.value  to:out_file			
		)				
	)	
		
-- Volume Integrator type	
fn Pbrt_VolumeIntegrator volType &out_file = (
		
	format "\n" 	to:out_file
		
	-- get settings
	if volType == 1 then --"single"
		(
			format "VolumeIntegrator \"single\" \n" to:out_file
			format "\t \"float stepsize\" [%] \n" PBRT_Render.Render_Vol_Step_spn.value to:out_file
		)		
	else if volType ==  2 then --"emission"
		(
			format "VolumeIntegrator \"emission\" \n" to:out_file
			format "\t \"float stepsize\" [%] \n" PBRT_Render.Render_Vol_Step_spn.value to:out_file
		)
	)	

-- Render type	
fn Pbrt_Renderer renType &out_file = (	

	format "\n" 	to:out_file
		
	-- get settings
	if renType == 1 then -- Aggregate Test
		(
			format "Renderer \"aggregatetest\" \n" to:out_file
			format "\t \"integer niters\" [%] \n" PBRT_Render.Render_Rnd_Aggre_Niters_spn.value to:out_file
		)		
	else if renType == 2 then -- Create Probes
		(
			format "Renderer \"createprobes\" \n" to:out_file
			format "\t \"bool directlighting\" [\"%\"] \n" PBRT_Render.Render_Rnd_Probe_Direct_btn.state to:out_file
			format "\t \"bool indirectlighting\" [\"%\"] \n" PBRT_Render.Render_Rnd_Probe_Indirect_btn.state to:out_file	
			format "\t \"integer lmax\" [%] \n" PBRT_Render.Render_Rnd_Probe_Lmax_spn.value to:out_file
			format "\t \"integer indirectsamples\" [%] \n" PBRT_Render.Render_Rnd_Probe_IndSamp_spn.value to:out_file			
			format "\t \"float samplespacing\" [%] \n" PBRT_Render.Render_Rnd_Probe_SampSpace_spn.value to:out_file
			format "\t \"float time\" [%] \n" PBRT_Render.Render_Rnd_Probe_time_spn.value to:out_file			
		)
	else if renType == 3 then -- Metropolis
		(
			format "Renderer \"metropolis\" \n" to:out_file
			format "\t \"float largestepprobability\" [%] \n" PBRT_Render.Render_Rnd_Metro_Probab_spn.value to:out_file
			format "\t \"integer nsamples\" [%] \n" PBRT_Render.Render_Rnd_Metro_Nsamples_spn.value to:out_file
			format "\t \"integer samplesperpixel\" [%] \n" PBRT_Render.Render_Rnd_Metro_Samples_spn.value to:out_file
			format "\t \"integer bootstrapsamples\" [%] \n" PBRT_Render.Render_Rnd_Metro_BootStr_spn.value to:out_file
			format "\t \"integer directsamples\" [%] \n" PBRT_Render.Render_Rnd_Metro_DirectSamp.value to:out_file
			format "\t \"bool dodirectseparately\" [\"%\"] \n" PBRT_Render.Render_Rnd_Metro_DoDirect.state to:out_file
			format "\t \"integer maxconsecutiverejects\" [%] \n" PBRT_Render.Render_Rnd_Metro_MaxReje_spn.value to:out_file
			format "\t \"integer maxdepth\" [%] \n" PBRT_Render.Render_Rnd_Metro_MaxDepth.value to:out_file		
			format "\t \"bool indirectonly\" [\"%\"] \n" PBRT_Render.Render_Rnd_Metro_Indirect_btn.state to:out_file			
		)
	else if renType == 4 then -- Standard
		(
			format "Renderer \"sample\" \n" to:out_file
		)
	else if renType == 5 then -- Surface Points
		(
			format "Renderer \"surfacepoints\" \n" to:out_file
			format "\t \"float minsampledistance\" [%] \n" PBRT_Render.Render_Rnd_Surface_MinSam_spn.value to:out_file			
		)		
	)
	
-- Textures
fn Pbrt_Textures obj type &out_file =
	(		
		Local map --= obj.material.KdTexture
		Local ExportedTex = #()
		
		Local textureString = #("","")

		Local MapType

		if type ==	"KdTexture" then (
			MapType = "imagemap" 
			map = obj.material.KdTexture
		)
			else if type == "bumpMap" then (
				MapType = "scale"
				map = obj.material.bumpMap
			) else if type == "sigma" then (
					MapType = "sigma"
					map = obj.material.sigmaMap				
				)else if type == "KsMap" then (
					MapType = "KsMap"
					map = obj.material.KsTexture				
					)else if type == "RoughMap" then (
						MapType = "RoughMap"
						map = obj.material.RoughnessTexture				
					)else if type == "KrMap" then (
							MapType = "KrMap"
							map = obj.material.KrTexture				
						)else if type == "OpacityMap" then (
							MapType = "OpacityMap"
							map = obj.material.opacityTexture				
						)else if type == "LightMap" then (
							MapType = "LightMap"
							map = obj.material.LightCMap				
						)else if type == "AlphaMap" then (
							MapType = "AlphaMap"
							map = obj.material.alphaImage				
						)else if type == "ReflectMap" then (
							MapType = "ReflectMap"
							map = obj.material.reflectTexture				
						)else if type == "TransmitMap" then (
							MapType = "TransmitMap"
							map = obj.material.transmitTexture				
						)else if type == "uRoughMap" then (
							MapType = "uRoughMap"
							map = obj.material.uRoughnessTexture			
						)else if type == "vRoughMap" then (
							MapType = "vRoughMap"
							map = obj.material.vRoughnessTexture			
						)

						
		if (classof map == PBRT_BitMap) then
			(
	
				local wrap 			= #("repeat", "black", "clamp") 	
				local mapType	= #("uv", "spherical", "cylindrical", "planar")
				local mapTypeID	= map.uv
				local wrapID 		= map.wrap
--				local scl 			= map._scale
				local file = undefined
				
				-- check if the user wants to use the default texture path for each texture or one specified by the user.
				-- this is usfull when sending the file to someone else, so all the maps are in the same place.
				if PBRT_Tools.uncPath_chk.checked then (
						Local theMap = filenameFromPath map.PBRT_BitMap.filename
					--	file = "\"" + PBRT_Tools.uncPath_etxt.text + "\\\\" + theMap + "\""
						file = "\"" + PBRT_Tools.uncPath_etxt.text + "/" + theMap + "\""
					)
						else (
								file = if map.PBRT_BitMap != undefined then PbrtAddBack map.PBRT_BitMap.filename
								)

				if type == "bumpMap" or type == "sigma" or type == "AlphaMap" or type == "RoughMap"  or type == "uRoughMap" or type == "vRoughMap" then (
					format "Texture \"%\" \"float\" \"imagemap\"\n" map.name to:out_file
				)else (
				format "Texture \"%\" \"color\" \"imagemap\"\n" map.name to:out_file
				)
				format "\t \"string mapping\" \"%\"\n" mapType[mapTypeID] to:out_file			
				format "\t \"string wrap\" [\"%\"] \n" wrap[wrapID] to:out_file
				format "\t \"string filename\" [%]\n"  file to:out_file
				format "\t \"float gamma\" [%]\n" map.gamma to:out_file
				format "\t \"string mapping\" [\"%\"]\n" mapType[mapTypeID] to:out_file
				format "\t \"float uscale\" [%]\n" map.uscale to:out_file
				format "\t \"float vscale\" [%]\n" map.vscale to:out_file
				format "\t \"float udelta\" [%]\n" map.udelta to:out_file
				format "\t \"float vdelta\" [%]\n" map.vdelta to:out_file	
				format "\t \"float maxanisotropy\" [%]\n" map.anist to:out_file
				format "\t \"bool trilinear\" [\"%\"]\n" map.tri to:out_file			

				if type == "bumpMap" then (
					Local MapName = map.name + "_bump"
					format "Texture \"%\" \"float\" \"scale\" \"texture tex1\" \"%\" \"float tex2\" [%]  \n" MapName map.name obj.material.bumpMapAmount to:out_file
				)
				
				if type == "sigma" then (
					Local MapName = map.name + "_sigma"
					format "Texture \"%\" \"float\" \"scale\" \"texture tex1\" [\"%\"] \"float tex2\" [%]  \n" MapName map.name obj.material.sigma to:out_file
				)
				
				if type == "KsMap" then (
					Local MapName = map.name + "_KsMap"
					format "Texture \"%\" \"color\" \"scale\" \"texture tex1\" [\"%\"] \"color tex2\" [% % %]  \n" MapName map.name 1.0 1.0 1.0 to:out_file
				)
				
				if type == "KrMap" then (
					Local MapName = map.name + "_KrMap"
					format "Texture \"%\" \"color\" \"scale\" \"texture tex1\" [\"%\"] \"color tex2\" [% % %]  \n" MapName map.name 1.0 1.0 1.0 to:out_file
				)
				
				if type == "OpacityMap" then (
					Local MapName = map.name + "_OpacityMap"
					format "Texture \"%\" \"color\" \"scale\" \"texture tex1\" [\"%\"] \n" MapName map.name  to:out_file
				)
				
				if type == "RoughMap" then (
					Local MapName = map.name + "_RoughMap"
					format "Texture \"%\" \"float\" \"scale\" \"texture tex1\" \"%\" \"float tex2\" [%]  \n" MapName map.name obj.material.roughness to:out_file
				)
				
				if type == "ReflectMap" then (
					Local MapName = map.name + "_ReflectMap"
					format "Texture \"%\" \"color\" \"scale\" \"texture tex1\" [\"%\"] \"color tex2\" [% % %]  \n" MapName map.name 1.0 1.0 1.0 to:out_file
				)		
				
				if type == "TransmitMap" then (
					Local MapName = map.name + "_TransmitMap"
					format "Texture \"%\" \"color\" \"scale\" \"texture tex1\" [\"%\"] \"color tex2\" [% % %]  \n" MapName map.name 1.0 1.0 1.0 to:out_file
				)
				
				if type == "uRoughMap" then (
					Local MapName = map.name + "_uRoughMap"
					format "Texture \"%\" \"float\" \"scale\" \"texture tex1\" \"%\" \"float tex2\" [%]  \n" MapName map.name obj.material.uRoughness to:out_file
				)
				
				if type == "vRoughMap" then (
					Local MapName = map.name + "_vRoughMap"
					format "Texture \"%\" \"float\" \"scale\" \"texture tex1\" \"%\" \"float tex2\" [%]  \n" MapName map.name obj.material.vRoughness to:out_file
				)				

				
/*					
				if type == "AlphaMap" then (
--					Local MapName = map.name + "_AlphaMap"
--					format "Texture \"%\" \"float\" \"scale\" \"texture tex1\" \"%\" \"float tex2\" [%]  \n" MapName map.name obj.material.bumpMapAmount to:out_file
				)
			
				if type == "LightMap" then (
					Local MapName = map.name + "_LightMap"
					format "Texture \"%\" \"color\" \"imagemap\" \"texture tex1\" [\"%\"] \"float tex2\" [%]  \n" MapName map.name obj.material.LightCMap to:out_file
				)
*/
--			)
			)
		else if(classof map == PBRT_CheckboardMap) then (

				local aa = #("closedform", "supersample", "none")
				local aaTypeID = map.aamode
				local mapType	= #("uv", "spherical", "cylindrical", "planar")
				local mapTypeID	= map.uv	

				format "Texture \"%\" \"color\" \"checkerboard\"\n" map.name  to:out_file
				format "\t \"string mapping\" \"%\"\n" mapType[mapTypeID] to:out_file
				format "\t \"color tex1\" [%] \n" (PbrtColor map.color1)  to:out_file
				format "\t \"color tex2\" [%] \n" (PbrtColor map.color2)  to:out_file
				format "\t \"string aamode\" [\"%\"]\n" aa[aaTypeID] to:out_file
				format "\t \"float uscale\" [%]\n" map.uscale to:out_file
				format "\t \"float vscale\" [%]\n" map.vscale to:out_file
				format "\t \"float udelta\" [%]\n" map.udelta to:out_file
				format "\t \"float vdelta\" [%]\n" map.vdelta to:out_file			


				if type == "bumpMap" then (
					Local MapName = uniquename (map.name + "_bump" )
					format "Texture \"%\" \"scale\" \"texture tex1\" \"%\" \"float tex2\" [%]  \n" MapName map.name obj.material.bumpMapAmount to:out_file
				)			
						
		)
		else if (classof map == PBRT_ConstantMap) then (

				format "Texture \"%\" \"color\" \"scale\"\n" map.name  to:out_file
				format "\t \"color tex1\" [%] \n" (PbrtColor map.color)  to:out_file
				format "\t \"color tex2\" [%] \n" (PbrtColor obj.material.kd)  to:out_file									
		)		
	)	
	
-- Materials		
fn Pbrt_Material obj  &out_file = (

		Local mat = obj.material
	
		-- first check if the user wants to use an overrride material
		if PBRT_Render.useGlobalMat_chk.checked and GlobalMaterial != undefined then mat = GlobalMaterial
		
		-- if it dosent have a material create one.
		if mat == undefined then (
				theMat = PBRT_Matte ()
				theMat.kd = (color 160 160 160)
				mat = theMat
			)

		-- check if the material hasnt been processed
		if not (finditem doneMats mat) != 0 and classof mat != PBRT_Light  then (
				append doneMats mat	-- add it to the array so we dont process it again	
						
				-- write object name
				format "\n### Object: % ### \n" Obj.name  to:out_file
					
		-- ************     Matte Material     ***********************
				if (classof mat == PBRT_Matte) then (
			
					local Kd = PbrtColor mat.Kd
	
					Local hasDiffuseMap = false
					Local hasSigmaMap	= false					
					Local hasBumpMap	= false
					Local hasAlphaMap	= false
--							Local hasLightMap	= false
					
					------- *** Create the Textures ***
					-- Check if we want to create the textures or not
					if PBRT_System.PBRT_DisableMaps_chk.checked != true then (
						
							-- Diffuse Texture
							if mat.KdTextureEnable == true and mat.KdTexture != undefined then (
								hasDiffuseMap = true
								
								-- first check if the texture hasnt been created jet
								if (finditem TextArray mat.KdTexture) == 0 then (
									Pbrt_Textures obj "KdTexture" out_file	-- create the texture
									append TextArray mat.KdTexture
									)
							)
							-- Sigma Texture
							if mat.sigmaMapEnabled == true and mat.sigmaMap != undefined then (
								hasSigmaMap = true
								-- first check if the texture hasnt been created jet
								if (finditem TextArray mat.sigmaMap) == 0 then (
									Pbrt_Textures obj "sigma" out_file -- create the texture
									append TextArray mat.sigmaMap
									)
							)					
							-- Bump Map Texture
							if mat.bumpMapEnable == true and mat.bumpMap != undefined then (
								hasBumpMap	= true
								
								-- first check if the texture hasnt been created jet
								if (finditem TextArray mat.bumpMap) == 0 then (
										Pbrt_Textures obj "bumpMap"	out_file -- create the texture
										append TextArray mat.bumpMap
									)
							)
							-- Alpha Map Texture
							if mat.alphaEnable == true and mat.alphaImage != undefined then (
								hasAlphaMap	= true
								
								-- first check if the texture hasnt been created jet
								if (finditem TextArray mat.alphaImage) == 0 then (
										Pbrt_Textures obj "AlphaMap"	out_file -- create the texture
										append TextArray mat.alphaImage
									)
							)									
																
						)
					
						
					-- *** Create The Material *** --

					format "MakeNamedMaterial \"%\" \n" mat.name  to:out_file
					format "\t \"string type\" [\"matte\"] \n" to:out_file		
					
					
					-- if it dosent ha a map or its not beenuse then use the color
					if mat.KdTextureEnable == false  or mat.KdTexture == undefined  then (
					format "\t \"color Kd\" [%]\n"  kd to:out_file
					)
					
					-- if the maps are diable use the diffuse color only
					if PBRT_System.PBRT_DisableMaps_chk.checked == true then format "\t \"color Kd\" [%]\n"  kd to:out_file
					
					if mat.sigma > 0.0 and mat.sigmaMapEnabled == false then format "\t \"float sigma\" [%]\n"  mat.sigma to:out_file
				
					if hasDiffuseMap then format "\t \"texture Kd\" [\"%\"]\n"  mat.KdTexture.name to:out_file
					
					if hasSigmaMap then (
						Local MapSigName = mat.sigmaMap.name + "_sigma"						
						format "\t \"texture sigma\" [\"%\"]\n"  MapSigName to:out_file
					)

					if hasBumpMap then (
						Local BumpMapName = mat.bumpMap.name + "_bump"						
						format "\t \"texture bumpmap\" [\"%\"]\n"  BumpMapName to:out_file							
					)	
						
				)				
		-- ************     Plastic Material     ***********************
				else if (classof mat == PBRT_Plastic) then (
						
					local Kd = PbrtColor mat.Kd
					local Ks = PbrtColor mat.Ks
					
					Local hasDiffuseMap = false
					Local hasBumpMap	= false
					Local hasKsMap		= false
					Local hasRoughMap	= false
					Local hasAlphaMap	= false 
					
					
					if PBRT_System.PBRT_DisableMaps_chk.checked != true then (
						------- Create the Textures
						if mat.KdTextureEnable == true and mat.KdTexture != undefined then (
							hasDiffuseMap = true
							
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.KdTexture) == 0 then (
								Pbrt_Textures obj "KdTexture" out_file	-- create the texture
								append TextArray mat.KdTexture
								)
						)
						
						if mat.KsTextureEnable == true and mat.KsTexture != undefined then (
							hasKsMap = true
							
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.KsTexture) == 0 then (
									Pbrt_Textures obj "KsMap" out_file	-- create the texture
									append TextArray mat.KsTexture
								)
						)					
						
						if mat.bumpMapEnable == true and mat.bumpMap != undefined then (
							hasBumpMap	= true
							
						-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.bumpMap) == 0 then (
									Pbrt_Textures obj "bumpMap"	out_file -- create the texture
									append TextArray mat.bumpMap
								)
						)
							
						if mat.RoughnessTextureEnable == true and mat.RoughnessTexture != undefined then (
							hasRoughMap	= true
							
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.RoughnessTexture) == 0 then (
										Pbrt_Textures obj "RoughMap" out_file	-- create the texture
										append TextArray mat.RoughnessTexture
								)
						)	
						
						-- Alpha Map Texture
						if mat.alphaEnable == true and mat.alphaImage != undefined then ( --**AQUI
								hasAlphaMap	= true
								
								-- first check if the texture hasnt been created jet
								if (finditem TextArray mat.alphaImage) == 0 then (
										Pbrt_Textures obj "AlphaMap"	out_file -- create the texture
										append TextArray mat.alphaImage
									)
							)									
					)							
					
					
					-- *** Create The Material *** --
					format "MakeNamedMaterial \"%\" \n" mat.name  to:out_file
					format "\t \"string type\" [\"plastic\"] \n" to:out_file
					if mat.KsTextureEnable == false then format "\t \"color Ks\" [%]\n"  ks to:out_file
					if mat.RoughnessTextureEnable == false then format "\t \"float roughness\"  [%] \n"  mat.Roughness to:out_file

					-- if it dosent have a map or its not beenuse then use the color
					if mat.KdTextureEnable == false or mat.KdTexture == undefined  then (
					format "\t \"color Kd\" [%]\n"  kd to:out_file
					)
					
					-- if the maps are diable use the diffuse color only
					if PBRT_System.PBRT_DisableMaps_chk.checked == true then format "\t \"color Kd\" [%]\n"  kd to:out_file
					
					if hasDiffuseMap then format "\t \"texture Kd\" [\"%\"]\n"  mat.KdTexture.name to:out_file
						
					
						if hasKsMap then (
							Local MapKsName = mat.KsTexture.name + "_KsMap"						
							format "\t \"texture Ks\" [\"%\"]\n"  MapKsName to:out_file
						)
						if hasRoughMap then (
							Local RoughMapName = mat.RoughnessTexture.name + "_RoughMap"						
							format "\t \"texture roughness\" \"%\"\n"  RoughMapName to:out_file								
						)
						if hasBumpMap then (
							Local BumpMapName = mat.bumpMap.name + "_bump"						
							format "\t \"texture bumpmap\" [\"%\"]\n"  BumpMapName to:out_file							
						)				
								
				)							
		-- ************     Uber Material     ***********************
				else if (classof mat == PBRT_Uber ) then (
						
					local Kd = PbrtColor mat.Kd
					local Ks = PbrtColor mat.Ks
					local Kr = PbrtColor mat.Kr
					local Opacity = PbrtColor mat.opacity
					
					Local hasDiffuseMap 	= false	
					Local hasKsMap			= false
					Local hasKrMap			= false
					Local hasOpacityMap	= false
					Local hasRoughMap		= false
					Local hasBumpMap		= false					
					
					------- Create the Textures
					if PBRT_System.PBRT_DisableMaps_chk.checked != true then (
						if mat.KdTextureEnable == true and mat.KdTexture != undefined then (
							hasDiffuseMap = true
							
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.KdTexture) == 0 then (
								Pbrt_Textures obj "KdTexture" out_file	-- create the texture
								append TextArray mat.KdTexture
								)
						)
						
						if mat.KsTextureEnable == true and mat.KsTexture != undefined then (
							hasKsMap = true
							
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.KsTexture) == 0 then (
								Pbrt_Textures obj "KsMap" out_file	-- create the texture
								append TextArray mat.KsTexture
								)
						)					
						
						if mat.KrTextureEnable == true and mat.KrTexture != undefined then (
							hasKrMap	= true
																
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.KrTexture) == 0 then (
								Pbrt_Textures obj "KrMap" out_file	-- create the texture
								append TextArray mat.KrTexture
								)
						)							
						
						if mat.opacityTextureEnable == true and mat.opacityTexture != undefined then (
							hasOpacityMap	= true
																									
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.opacityTexture) == 0 then (
								Pbrt_Textures obj "OpacityMap" out_file	-- create the texture
								append TextArray mat.opacityTexture
								)
						)
						
						if mat.bumpMapEnable == true and mat.bumpMap != undefined then (
							hasBumpMap	= true
																																		
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.bumpMap) == 0 then (
								Pbrt_Textures obj "bumpMap" out_file	-- create the texture
								append TextArray mat.bumpMap
								)
						)
							
						if mat.RoughnessTextureEnable == true and mat.RoughnessTexture != undefined then (
							hasRoughMap	= true
							
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.RoughnessTexture) == 0 then (
								Pbrt_Textures obj "RoughMap" out_file	-- create the texture
								append TextArray mat.RoughnessTexture
								)
						)	
					)							
					
					
					-- *** Create The Material *** --
					format "MakeNamedMaterial \"%\" \n" mat.name  to:out_file
					format "\t \"string type\" [\"uber\"] \n" to:out_file
					
					if mat.KsTextureEnable == false then format "\t \"color Ks\"  [%] \n"  ks to:out_file
					if mat.KrTextureEnable == false then format "\t \"color Kr\"  [%] \n"  kr to:out_file
					if mat.opacityTextureEnable == false then format "\t \"color opacity\"  [%] \n"  opacity to:out_file
					if mat.RoughnessTextureEnable == false then format "\t \"float roughness\"  [%] \n"  mat.Roughness to:out_file

					-- if it dosent has a map or its not been use then use the color
					if mat.KdTextureEnable == false or mat.KdTexture == undefined  then (
					format "\t \"color Kd\" [%]\n"  kd to:out_file
					)
					
					-- if the maps are diable use the diffuse color only
					if PBRT_System.PBRT_DisableMaps_chk.checked == true then format "\t \"color Kd\" [%]\n"  kd to:out_file
					
						if hasDiffuseMap then format "\t \"texture Kd\" \"%\" \n"  mat.KdTexture.name to:out_file
							
						if hasKsMap then (
							Local MapKsName = mat.KsTexture.name + "_KsMap"						
							format "\t \"texture Ks\" \"%\"\n"  MapKsName to:out_file
						)
						if hasKrMap then (
							Local MapKrName = mat.KrTexture.name + "_KrMap"						
							format "\t \"texture Kr\" \"%\"\n"  MapKrName to:out_file
						)								
						if hasRoughMap then (
							Local RoughMapName = mat.RoughnessTexture.name + "_RoughMap"						
							format "\t \"texture roughness\" \"%\"\n"  RoughMapName to:out_file							
						)
						if hasOpacityMap then (
							Local OpacityMapName = mat.opacityTexture.name + "_OpacityMap"						
							format "\t \"texture opacity\" \"%\"\n"  OpacityMapName to:out_file							
						)								
						if hasBumpMap then (
							Local BumpMapName = mat.bumpMap.name + "_bump"						
							format "\t \"texture bumpmap\" \"%\"\n"  BumpMapName to:out_file							
						)				
								
				)							
		-- ************     Glass Material     ***********************
				else if (classof mat == PBRT_Glass) then (
					
					local Kr = PbrtColor mat.Kr
					local Kt = PbrtColor mat.Kt

					format "MakeNamedMaterial \"%\" \n" mat.name  to:out_file
					format "\t \"string type\" [\"glass\"] \n" to:out_file		

					
					Local hasBumpMap	= false
			
					------- Create the Textures				
					if PBRT_System.PBRT_DisableMaps_chk.checked != true then (
					if mat.bumpMapEnable == true and mat.bumpMap != undefined then (
								hasBumpMap	= true
						
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.bumpMap) == 0 then (
								Pbrt_Textures obj "bumpMap" out_file	-- create the texture
								append TextArray mat.bumpMap
								)
							)				
						)
					
					-- *** Create The Material *** --
					
						format "\t \"color Kr\" [%] \"color Kt\" [%] \"float index\" [%] \n"  Kr Kt mat.ior to:out_file	
					
						if hasBumpMap then (
							Local BumpMapName = mat.bumpMap.name + "_bump"						
							format "\t \"texture bumpmap\" [\"%\"]\n"  BumpMapName to:out_file							
						)				
										
				)				
		-- ************    Translucent Material     ***********************
				else if (classof mat == PBRT_Translucent) then (
						
					local Kd = PbrtColor mat.Kd
					local Ks = PbrtColor mat.Ks
					
					Local hasDiffuseMap = false
					Local hasBumpMap	= false
					Local hasReflectMap	= false
					Local hasTransmitMap	= false
					Local hasKsMap		= false
					Local hasRoughMap	= false
					Local hasAlphaMap	= false 
					Local hasAlphaMap	= false 
					
					
				if PBRT_System.PBRT_DisableMaps_chk.checked != true then (
						------- Create the Textures
						if mat.KdTextureEnable == true and mat.KdTexture != undefined then (
							hasDiffuseMap = true
							
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.KdTexture) == 0 then (
								Pbrt_Textures obj "KdTexture" out_file	-- create the texture
								append TextArray mat.KdTexture
								)
						)
						
						if mat.KsTextureEnable == true and mat.KsTexture != undefined then (
							hasKsMap = true
							
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.KsTexture) == 0 then (
									Pbrt_Textures obj "KsMap" out_file	-- create the texture
									append TextArray mat.KsTexture
								)
						)					
						
						if mat.bumpMapEnable == true and mat.bumpMap != undefined then (
							hasBumpMap	= true
							
						-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.bumpMap) == 0 then (
									Pbrt_Textures obj "bumpMap"	out_file -- create the texture
									append TextArray mat.bumpMap
								)
						)
							
						if mat.RoughnessTextureEnable == true and mat.RoughnessTexture != undefined then (
							hasRoughMap	= true
							
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.RoughnessTexture) == 0 then (
										Pbrt_Textures obj "RoughMap" out_file	-- create the texture
										append TextArray mat.RoughnessTexture
								)
						)	
							
						if mat.reflectTextureEnable == true and mat.reflectTexture != undefined then (
							hasReflectMap	= true
							
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.reflectTexture) == 0 then (
										Pbrt_Textures obj "ReflectMap" out_file	-- create the texture
										append TextArray mat.reflectTexture
								)
						)							
							
						if mat.transmitTextureEnable == true and mat.transmitTexture != undefined then (
							hasTransmitMap	= true
							
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.transmitTexture) == 0 then (
										Pbrt_Textures obj "TransmitMap" out_file	-- create the texture
										append TextArray mat.transmitTexture
								)
						)	
						
						-- Alpha Map Texture
						if mat.alphaEnable == true and mat.alphaImage != undefined then ( --**AQUI
								hasAlphaMap	= true
								
								-- first check if the texture hasnt been created jet
								if (finditem TextArray mat.alphaImage) == 0 then (
										Pbrt_Textures obj "AlphaMap"	out_file -- create the texture
										append TextArray mat.alphaImage
									)
							)
							
					)							
					
					
					-- *** Create The Material *** --
					format "MakeNamedMaterial \"%\" \n" mat.name  to:out_file
					format "\t \"string type\" [\"translucent\"] \n" to:out_file
					if mat.KsTextureEnable == false then format "\t \"color Ks\" [%]\n"  ks to:out_file
					if mat.RoughnessTextureEnable == false then format "\t \"float roughness\"  [%] \n"  mat.Roughness to:out_file
					if mat.reflectTextureEnable == false then format "\t \"color reflect\" [% % %]\n"  mat.reflect mat.reflect mat.reflect to:out_file
					if mat.transmitTextureEnable == false then format "\t \"color transmit\" [% % %]\n"  mat.transmit mat.transmit mat.transmit to:out_file
						
					-- if it dosent have a map or its not beenuse then use the color
					if mat.KdTextureEnable == false or mat.KdTexture == undefined  then (
					format "\t \"color Kd\" [%]\n"  kd to:out_file
					)
					
					-- if the maps are diable use the diffuse color only
					if PBRT_System.PBRT_DisableMaps_chk.checked == true then format "\t \"color Kd\" [%]\n"  kd to:out_file
					
					if hasDiffuseMap then format "\t \"texture Kd\" [\"%\"]\n"  mat.KdTexture.name to:out_file
						
					
						if hasKsMap then (
							Local MapKsName = mat.KsTexture.name + "_KsMap"						
							format "\t \"texture Ks\" [\"%\"]\n"  MapKsName to:out_file
						)
						if hasRoughMap then (
							Local RoughMapName = mat.RoughnessTexture.name + "_RoughMap"						
							format "\t \"texture roughness\" \"%\"\n"  RoughMapName to:out_file								
						)
						if hasBumpMap then (
							Local BumpMapName = mat.bumpMap.name + "_bump"						
							format "\t \"texture bumpmap\" [\"%\"]\n"  BumpMapName to:out_file							
						)	
						if hasReflectMap then (
							Local ReflectName = mat.reflectTexture.name + "_ReflectMap"						
							format "\t \"texture reflect\" [\"%\"]\n"  ReflectName to:out_file
						)
						if hasTransmitMap then (
							Local TransmitName = mat.transmitTexture.name + "_TransmitMap"						
							format "\t \"texture transmit\" [\"%\"]\n"  TransmitName to:out_file
						)						
						
				)	
		-- ************    Substrate Material     ***********************
				else if (classof mat == PBRT_Substrate) then (
					local Kd = PbrtColor mat.Kd
					local Ks = PbrtColor mat.Ks
					
					Local hasDiffuseMap = false
					Local hasKsMap		= false
					Local hasURoughMap	= false
					Local hasVRoughMap	= false
					Local hasBumpMap	= false
					
					
					if PBRT_System.PBRT_DisableMaps_chk.checked != true then (
						------- Create the Textures
						if mat.KdTextureEnable == true and mat.KdTexture != undefined then (
							hasDiffuseMap = true
							
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.KdTexture) == 0 then (
								Pbrt_Textures obj "KdTexture" out_file	-- create the texture
								append TextArray mat.KdTexture
								)
						)
						
						if mat.KsTextureEnable == true and mat.KsTexture != undefined then (
							hasKsMap = true
							
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.KsTexture) == 0 then (
									Pbrt_Textures obj "KsMap" out_file	-- create the texture
									append TextArray mat.KsTexture
								)
						)					
						
						if mat.bumpMapEnable == true and mat.bumpMap != undefined then (
							hasBumpMap	= true
							
						-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.bumpMap) == 0 then (
									Pbrt_Textures obj "bumpMap"	out_file -- create the texture
									append TextArray mat.bumpMap
								)
						)
							
						if mat.uRoughnessTextureEnable == true and mat.uRoughnessTexture != undefined then (
							hasURoughMap	= true
							
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.uRoughnessTexture) == 0 then (
										Pbrt_Textures obj "uRoughMap" out_file	-- create the texture
										append TextArray mat.uRoughnessTexture
								)
						)	
							
						if mat.vRoughnessTextureEnable == true and mat.vRoughnessTexture != undefined then (
							hasVRoughMap	= true
							
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.vRoughnessTexture) == 0 then (
										Pbrt_Textures obj "vRoughMap" out_file	-- create the texture
										append TextArray mat.vRoughnessTexture
								)
						)							
								
					)							
					
					
					-- *** Create The Material *** --
					format "MakeNamedMaterial \"%\" \n" mat.name  to:out_file
					format "\t \"string type\" [\"substrate\"] \n" to:out_file
					if mat.KsTextureEnable == false then format "\t \"color Ks\" [%]\n"  ks to:out_file
						
					if mat.uRoughnessTextureEnable == false then format "\t \"float uroughness\"  [%] \n"  mat.uRoughness to:out_file
					if mat.vRoughnessTextureEnable == false then format "\t \"float vroughness\"  [%] \n"  mat.vRoughness to:out_file

					-- if it dosent have a map or its not beenuse then use the color
					if mat.KdTextureEnable == false or mat.KdTexture == undefined  then (
					format "\t \"color Kd\" [%]\n"  kd to:out_file
					)
					
					-- if the maps are diable use the diffuse color only
					if PBRT_System.PBRT_DisableMaps_chk.checked == true then format "\t \"color Kd\" [%]\n"  kd to:out_file
					
					if hasDiffuseMap then format "\t \"texture Kd\" [\"%\"]\n"  mat.KdTexture.name to:out_file
						
					
						if hasKsMap then (
							Local MapKsName = mat.KsTexture.name + "_KsMap"						
							format "\t \"texture Ks\" [\"%\"]\n"  MapKsName to:out_file
						)
						if hasURoughMap then (
							Local uRoughMapName = mat.uRoughnessTexture.name + "_uRoughMap"						
							format "\t \"texture uroughness\" \"%\"\n"  uRoughMapName to:out_file								
						)
						if hasVRoughMap then (
							Local vRoughMapName = mat.vRoughnessTexture.name + "_vRoughMap"						
							format "\t \"texture vroughness\" \"%\"\n"  vRoughMapName to:out_file								
						)						
						if hasBumpMap then (
							Local BumpMapName = mat.bumpMap.name + "_bump"						
							format "\t \"texture bumpmap\" [\"%\"]\n"  BumpMapName to:out_file							
						)					
				)
		-- ************    Metal Material     ***********************
				else if (classof mat == PBRT_Metal) then (
						
					Local hasRoughMap	= false
					Local hasBumpMap	= false
					
					
				if PBRT_System.PBRT_DisableMaps_chk.checked != true then (
						------- Create the Textures				
						
						if mat.bumpMapEnable == true and mat.bumpMap != undefined then (
							hasBumpMap	= true
							
						-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.bumpMap) == 0 then (
									Pbrt_Textures obj "bumpMap"	out_file -- create the texture
									append TextArray mat.bumpMap
								)
						)
							
						if mat.RoughnessTextureEnable == true and mat.RoughnessTexture != undefined then (
							hasRoughMap	= true
							
							-- first check if the texture hasnt been created jet
							if (finditem TextArray mat.RoughnessTexture) == 0 then (
										Pbrt_Textures obj "RoughMap" out_file	-- create the texture
										append TextArray mat.RoughnessTexture
								)
						)										
					)							
					
					
					-- *** Create The Material *** --
					format "MakeNamedMaterial \"%\" \n" mat.name  to:out_file
					format "\t \"string type\" [\"metal\"] \n" to:out_file
					
					if mat.RoughnessTextureEnable == false then format "\t \"float roughness\"  [%] \n"  mat.Roughness to:out_file

					
					-- check if the user wants to use the default texture path for each texture or one specified by the user.
					-- this is usfull when sending the file to someone else, so all the maps are in the same place.
					if PBRT_Tools.UNCPath_chk.checked then (
								Local theEta = filenameFromPath mat.eta
								Local theK = filenameFromPath mat.k
					
						--		fileEta 	= PBRT_Tools.uncPath_etxt.text + "\\\\" + theEta
						--		fileK 		= PBRT_Tools.uncPath_etxt.text + "\\\\" + theK						
								fileEta 	= PBRT_Tools.uncPath_etxt.text + "/" + theEta
								fileK 		= PBRT_Tools.uncPath_etxt.text + "/" + theK
						
								format "\t \"spectrum eta\" \"%\" \n" fileEta to:out_file
								format "\t \"spectrum k\" \"%\" \n" 	fileK to:out_file
							)else (										
									format "\t \"spectrum eta\" % \n" (PbrtAddBack mat.eta) to:out_file
									format "\t \"spectrum k\" % \n" (PbrtAddBack mat.k) to:out_file
									)


						if hasRoughMap then (
							Local RoughMapName = mat.RoughnessTexture.name + "_RoughMap"						
							format "\t \"texture roughness\" \"%\"\n"  RoughMapName to:out_file								
						)
						if hasBumpMap then (
							Local BumpMapName = mat.bumpMap.name + "_bump"						
							format "\t \"texture bumpmap\" [\"%\"]\n"  BumpMapName to:out_file							
						)				
									
				)				
		-- ********** Any other Material ***************
				else
				(
					format "MakeNamedMaterial \"%\" \n" mat.name  to:out_file
					format "\t \"string type\" [\"matte\"] \n" to:out_file							
					format "\t \"color Kd\" [.75 .75 .75]\n"  kd to:out_file
				)
				
				format "\n### End Object: % ### \n" Obj.name  to:out_file
			)
				
)

-- Geometry Array attributes
fn Pbrt_GetMeshArrays Obj = (
	
			-- Array to hold the mesh arrays
	
			-- MeshArrays [1] Faces Array		
			-- MeshArrays [2] Vertex Array
			-- MeshArrays [3] Texture Vertex Array
			-- MeshArrays [4] Normal Vertex Array
			Local MeshArrays = #(#(), #(), #(), #())

			num_faces = Obj.numfaces
				
			for i = 1 to num_faces do
					(
						face = getFace Obj i
						
						v1 = (MeshArrays[2].count + 1)
						v2 = (MeshArrays[2].count + 2)
						v3 = (MeshArrays[2].count + 3)
						
						append MeshArrays[1] [v1,v2,v3]
	/*					
						v1 = getVert Obj face.x
						v2 = getVert Obj face.y
						v3 = getVert Obj face.z
	*/
						v1 = coordsys world getvert Obj face.x	
						v2 = coordsys world getvert Obj face.y
						v3 = coordsys world getvert Obj face.z
						
						append MeshArrays[2] v1
						append MeshArrays[2] v2
						append MeshArrays[2] v3
				 
						v1 = (getnormal Obj face.x) --* theInvTM 
						v2 = (getnormal Obj face.y) --* theInvTM 
						v3 = (getnormal Obj face.z) --* theInvTM 
					
						append MeshArrays[4] v1
						append MeshArrays[4] v2
						append MeshArrays[4] v3	
						
						if Obj.numtverts != 0 then ( -- check if it has uvw's
							tvface = getTVFace Obj i

							v1 = getTVert Obj tvface.x
							v2 = getTVert Obj tvface.y
							v3 = getTVert Obj tvface.z

							append MeshArrays[3] v1
							append MeshArrays[3] v2	
							append MeshArrays[3] v3	
							)
					)
			
		return MeshArrays
)

-- Instance Geometry
fn Pbrt_GeometryInstance  Obj &out_file = (
	
	-- first object

	format "\n" to:out_file   
	format "# Main Instance Object:%  \n" Obj[1].name to:out_file  
	format "\n" to:out_file  
	   
	format "ObjectBegin \"%\" \n" Obj[1].name to:out_file 
	format "\t AttributeBegin \n"  to:out_file 		

	oriTM = Obj[1].transform
	--tm = Obj[1].transform * inverse Obj[1].transform			
	tm = Obj[1].transform * inverse Obj[1].transform
	
	format "\n" to:out_file   
	format "\t Transform \n\t [\n" to:out_file 
	format "\t\t % % % 0 \n" tm[1][1] tm[1][2] tm[1][3]  to:out_file 
	format "\t\t % % % 0 \n" tm[2][1] tm[2][2] tm[2][3]   to:out_file 
	format "\t\t % % % 0 \n" tm[3][1] tm[3][2] tm[3][3]   to:out_file 
	format "\t\t % % % 1 \n" tm[4][1] tm[4][2] tm[4][3]   to:out_file 
	format "\t ]\n" to:out_file 
	

							
	format "\t NamedMaterial \"%\" \n" Obj[1].material.name to:out_file 
	format "\t Shape \"trianglemesh\" \n" to:out_file 
	format "\t \"integer indices\" \n" to:out_file 
	format "\t [\n" to:out_file 
	

	Local MeshArrays = #(#(), #(), #(), #())

	num_faces = Obj[1].numfaces
				
	for i = 1 to num_faces do
			(
				face = getFace Obj[1] i
				
				v1 = (MeshArrays[2].count + 1)
				v2 = (MeshArrays[2].count + 2)
				v3 = (MeshArrays[2].count + 3)
				
				append MeshArrays[1] [v1,v2,v3]

				v1 = coordsys local getvert Obj[1] face.x	
				v2 = coordsys local getvert Obj[1] face.y
				v3 = coordsys local getvert Obj[1] face.z
				
				append MeshArrays[2] v1
				append MeshArrays[2] v2
				append MeshArrays[2] v3
		 
				v1 = (getnormal Obj[1] face.x) --* theInvTM 
				v2 = (getnormal Obj[1] face.y) --* theInvTM 
				v3 = (getnormal Obj[1] face.z) --* theInvTM 
			
				append MeshArrays[4] v1
				append MeshArrays[4] v2
				append MeshArrays[4] v3	
				
				if Obj[1].numtverts != 0 then ( -- check if it has uvw's
					tvface = getTVFace Obj[1] i

					v1 = getTVert Obj[1] tvface.x
					v2 = getTVert Obj[1] tvface.y
					v3 = getTVert Obj[1] tvface.z

					append MeshArrays[3] v1
					append MeshArrays[3] v2	
					append MeshArrays[3] v3	
					)
			)

					-- Face Index
		for f in MeshArrays[1] do
					(
						format "\t% % %\n" (f.x as integer - 1) (f.y as integer - 1) (f.z as integer - 1) to:out_file
					)
					
		format "\t]\n" to:out_file 			
					
		-- Vertex Points	
		format "\t\"point P\"\n" to:out_file 					
		format "\t[\n" to:out_file 					
		for v in MeshArrays[2] do
					(			
						v.x = v.x -- Obj[1].transform.row4.x
						v.y = v.y -- Obj[1].transform.row4.y
						v.z = v.z -- Obj[1].transform.row4.z						
						format "\t\t% % %\n" v.x v.y v.z to:out_file
					)
		format "\t]\n"  to:out_file 
	
		-- just for matte materials for now.
		if classof Obj[1].material == PBRT_Matte or classof Obj[1].material == PBRT_Plastic or classof Obj[1].material == PBRT_Translucent then (
			
				if Obj[1].material.alphaAmount < 1.0 and Obj[1].material.alphaEnable == false then (
						format "\t \"float alpha\" [%]\n" Obj[1].material.alphaAmount  to:out_file	
					)
				if Obj[1].material.alphaEnable and Obj[1].material.alphaImage != undefined then (
						format "\t \"texture alpha\" \"%\"\n" Obj[1].material.alphaImage.name to:out_file	
					)									
					
			)
		
		-- **** Uv's information ****
		-- skip if there are no uv's
		if MeshArrays[3].count != 0 then (
			format "\t\"float uv\"\n"  to:out_file		
			format "\t[\n"  to:out_file
			for v in MeshArrays[3] do
			(
				vert = v
				format "\t\t% %\n" vert.x vert.y to:out_file
			)
			format "\t]\n"  to:out_file					
		)
		-- **** Normals ****
		format "\t\"normal N\"\n"  to:out_file	
		format "\t[\n"  to:out_file
		for v in MeshArrays[4] do
			(
				format "\t\t% % %\n" v.x v.y v.z to:out_file
			)
		format "\t]\n"  to:out_file
					
		format "\t AttributeEnd\n" to:out_file 
	   format "ObjectEnd \n" to:out_file 			
		
		-- rest of the objects	
					
		for o = 1 to Obj.count do (	
			
			TM = Obj[o].transform --* (inverse Obj[o].transform)
			
			format "\n" to:out_file   
			format "# Instanced Object:% \n" Obj[o].name to:out_file  
			format "AttributeBegin\n" to:out_file  
		
			format "\n" to:out_file   
			format "\t Transform \n\t [\n" to:out_file 
			format "\t\t % % % 0 \n" tm[1][1] tm[1][2] tm[1][3]  to:out_file 
			format "\t\t % % % 0 \n" tm[2][1] tm[2][2] tm[2][3]   to:out_file 
			format "\t\t % % % 0 \n" tm[3][1] tm[3][2] tm[3][3]   to:out_file 
			format "\t\t % % % 1 \n" tm[4][1] tm[4][2] tm[4][3]   to:out_file 
			format "\t ]\n" to:out_file 
/*	
			format "\t Transform \n\t [\n" to:out_file 
			format "\t\t % % % 0 \n" 1.0 0.0 0.0  to:out_file 
			format "\t\t % % % 0 \n" 0.0 1.0 0.0   to:out_file 
			format "\t\t % % % 0 \n" 0.0 0.0 1.0   to:out_file 
			format "\t\t % % % 1 \n" TM[4][1] TM[4][2] TM[4][3]   to:out_file 
			format "\t ]\n" to:out_file 

		
			trans = Obj[o].transform.translationpart -- oriTM.translationpart
			format "\t Translate  % % % \n" trans.x trans.y trans.z to:out_file 
		
			ang = Obj[o].transform.rotationpart.angle -- oriTM.rotationpart.angle
			axis = Obj[o].transform.rotationpart.axis * -1
				
--			format "\t Rotate  % % % % \n" ang axis.x axis.y axis.z to:out_file
			format "\t Rotate  % % % % \n" ang 0 0 1 to:out_file
	*/		
			format "\t ObjectInstance \"%\" \n" Obj[1].name to:out_file 
			format "AttributeEnd \n" to:out_file 			
		)
		
)

fn Pbrt_Geometry Obj  &out_file = (

					--- Object Geometry Output
				format "\n#**** Object: % **** \n"  Obj.name to:out_file
				format "AttributeBegin \n"  to:out_file
				format "\t NamedMaterial \"%\" \n" Obj.material.name  to:out_file

				-- i need to modify this, but for now just skip the spheres.
				if classof Obj != Sphere then  MeshArrays = Pbrt_GetMeshArrays Obj 

				-- Geometry Area Lights
				if classof Obj.material == PBRT_Light or \
					(classof Obj.material == PBRT_Matte and Obj.material.LightColorEnable) or \
						(classof Obj.material == PBRT_Uber and Obj.material.LightColorEnable) then (
							

							
					theColor = (PbrtColor Obj.material.LightColor)
					theColor2 = (PbrtColor (Obj.material.LightColor * 5))

					format "\t AreaLightSource \"area\" \"color L\" [%]\n" theColor to:out_file
					format "\t \"color scale\" [% % %]\n" Obj.material.lightMult Obj.material.lightMult Obj.material.lightMult to:out_file	
					format "\t \"integer nsamples\" [%] \n" Obj.material.LightSamples to:out_file
				
					Local theName = uniquename (Obj.name + "_Mat")
					
					format "\t Texture \"%\" \"color\" \"scale\"\n" theName to:out_file
					format "\t \"color tex1\" [%] \n" theColor2 to:out_file
--					format "\t \"color tex2\" [%] \n" theColor  to:out_file	
					format "\n" 	to:out_file						
					if classof Obj.material == PBRT_Matte then (
					if Obj.material.opacity == 255 then format "\n Material \"matte\" \"color Kd\" [%] \"texture Kd\" [\"%\"]\n" theColor theName to:out_file
					)
					
				)

				
				if classof Obj == Sphere then (
					
						format "\t Transform \n\t [\n" to:out_file 
						format "\t\t % % % 0 \n" Obj.transform[1][1] Obj.transform[1][2] Obj.transform[1][3]  to:out_file 
						format "\t\t % % % 0 \n" Obj.transform[2][1] Obj.transform[2][2] Obj.transform[2][3]   to:out_file 
						format "\t\t % % % 0 \n" Obj.transform[3][1] Obj.transform[3][2] Obj.transform[3][3]   to:out_file 
						format "\t\t % % % 1 \n" Obj.transform[4][1] Obj.transform[4][2] Obj.transform[4][3]   to:out_file 
						format "\t ]\n" to:out_file 
					
						format "\nShape \"sphere\" \"float radius\" [%]\n"  obj.radius to:out_file				
												
						format "\nAttributeEnd \n" to:out_file
						format "\n#**** End Object: % **** \n"  Obj.name to:out_file	
				)else (
						-- **** Triangle Indiex *****
						format "\nShape \"trianglemesh\" \"integer indices\"\n"  to:out_file	
						format "[\n"  to:out_file					
						for f in MeshArrays[1] do
							(
								format "\t% % %\n" (f.x as integer - 1) (f.y as integer - 1) (f.z as integer - 1) to:out_file
							)					
						format "]\n"  to:out_file
							
						-- **** Points ****
						format "\"point P\" \n[ \n"  to:out_file					
						for v in MeshArrays[2] do
							(
								format "\t% % %\n" v.x v.y v.z to:out_file
							)					
						format "]\n"  to:out_file

						-- Alpha Images or amount
							
						-- just for matte materials for now.
						if classof Obj.material == PBRT_Matte or classof Obj.material == PBRT_Plastic or classof Obj.material == PBRT_Translucent then (
							
								if Obj.material.alphaAmount < 1.0 and Obj.material.alphaEnable == false then (
										format "\t \"float alpha\" [%]\n" Obj.material.alphaAmount  to:out_file	
									)
								if Obj.material.alphaEnable and Obj.material.alphaImage != undefined then (
										format "\t \"texture alpha\" \"%\"\n" Obj.material.alphaImage.name to:out_file	
									)																		
							)
							
						if classof Obj.material ==  PBRT_Uber or classof Obj.material == PBRT_Light then (
								if Obj.material.emitterInv == true then format "\t \"float alpha\" [0]\n" to:out_file	
							
							)							
						----
							
						-- **** Uv's information ****
	--					if Obj.numtverts != 0 then ( -- check if it has uvw's
						if MeshArrays[3].count != 0 then (
							format "\"float uv\"  \n[\n"  to:out_file		
							for v in MeshArrays[3] do
							(
								vert = v -- Obj.transform.row4
								format "\t% %\n" vert.x vert.y to:out_file
							)
							format "]\n"  to:out_file
						)
/*	*/						
						-- **** Normals ****
						format "\"normal N\" \n[\n"  to:out_file	
						for v in MeshArrays[4] do
							(
								format "\t% % %\n" v.x v.y v.z to:out_file
							)
						format "]\n"  to:out_file
						
						format "\nAttributeEnd \n" to:out_file
						format "\n#**** End Object: % **** \n"  Obj.name to:out_file		
						)	

)
	
fn Pbrt_Lights Lt &out_file = (

	-- Point Lights
	if classof Lt == Omnilight then (
			ColMult = Lt.color * Lt.multiplier					
			theColor = (PbrtColor ColMult)

			format "# Point Light-> %\n"  Lt.name to:out_file		
			format "TransformBegin \n"  to:out_file
			format "\t Translate % % %\n" Lt.pos.x Lt.pos.y Lt.pos.z to:out_file
			format "\t LightSource \"point\" \n"  to:out_file
			format "\t \"point from\" [0 0 0] \n"  to:out_file 
			format "\t \"color I\" [%] \n" theColor to:out_file
			format "TransformEnd \n"	 Lt.name to:out_file	
			format "# End Point Light-> %\n\n"	 Lt.name to:out_file	
	)
	
	-- Spot Lights
	else if classof Lt == targetSpot  then (
		
/*			
			local LtNormal = Lt 
			local worldUpVector = [0,0,1]
			local rightVector = normalize (cross worldUpVector faceNormal)
			local upVector = normalize ( cross rightVector faceNormal )
			local theMatrix = matrix3 rightVector upVector faceNormal [0,0,0]
			local theDir = theMatrix.row3 * -1
*/		
			
			ColMult = Lt.color * Lt.multiplier					
			theColor = (PbrtColor ColMult)
		
			theDir = normalize (Lt.target.pos - Lt.pos)
		--	theDir = getDirection Lt
			format "# Spot Light-> %\n"  Lt.name to:out_file	
			format "TransformBegin		\n"  to:out_file
 				
			format "\t Translate % % %\n" Lt.pos.x Lt.pos.y Lt.pos.z to:out_file
			format "\t LightSource \"spot\" \n"  to:out_file
			format "\t \"point from\" [0 0 0] \n"  to:out_file
			format "\t \"point to\" [% % %] \n"  theDir.x theDir.y theDir.z to:out_file
			format "\t \"color I\" [%] \n" theColor to:out_file
			format "\t \"float coneangle\" [%] \n" Lt.hotspot to:out_file
			format "\t \"float conedeltaangle\" [%] \n" (Lt.falloff - Lt.hotspot) to:out_file
			format "TransformEnd \n"	 to:out_file	
			format "# End Spot Light-> %\n\n"	 Lt.name to:out_file	
	)	
	else if classof Lt == TargetDirectionallight then (
			ColMult = Lt.color * Lt.multiplier					
			theColor = (PbrtColor ColMult)
		
		
			theDir = normalize (Lt.target.pos - Lt.pos)
			format "# Direct Light-> %\n"  Lt.name to:out_file
			format "TransformBegin	\n"  to:out_file
 				
			format "\t Translate % % %\n" Lt.pos.x Lt.pos.y Lt.pos.z to:out_file
			format "\t LightSource \"distant\" \n"  to:out_file
			format "\t \"point from\" [0 0 0] \n"  to:out_file
			format "\t \"point to\" [% % %] \n"  theDir.x theDir.y theDir.z to:out_file
			format "\t \"color L\" [%] \n" theColor to:out_file
			format "TransformEnd \n" to:out_file		
			format "# End Direct Light-> %\n\n"	 Lt.name to:out_file	
	)
	
)

fn Pbrt_Accelerator Type &out_file= (
	
	if Type == 1 then -- None
		(
			--format "\n Accelerator \"\" % \n" PBRT_System. to:out_file
		)
	if Type == 2 then -- Grid
		(
			format "\n Accelerator \"grid\" \n"  to:out_file
			format "\t bool refineimmediately [\"%\"] \n" PBRT_System.Acce_Grid_Refine_btn.state to:out_file
			format "\n"   to:out_file
		)
	if Type == 3 then -- KdTree
		(
			format "\n Accelerator \"kdtreeaccel\" \n"  to:out_file
			format "\t \"integer intersectcost\" [%] \n" PBRT_System.Acce_KdTree_Intersect_spn.value to:out_file
			format "\t \"integer traversalcost\" [%] \n" PBRT_System.Acce_KdTree_Traversal_spn.value to:out_file
			format "\t \"float emptybonus\" [%] \n" PBRT_System.Acce_KdTree_EmptyBonus_Spn.value to:out_file
			format "\t \"integer maxprims\" [%] \n" PBRT_System.Acce_KdTree_MaxPrimes_spn.value to:out_file
			format "\t \"integer maxdepth\" [%] \n" PBRT_System.Acce_KdTree_MaxDepth_spn.value to:out_file
			format "\n"   to:out_file
		)
	if Type == 4 then -- bvh
		(
			format "\n Accelerator \"bvh\" \n"  to:out_file
		)		
)

-- Function to prepare the instance objects. for some reason they need to have a turn to mesh modifier??
fn prepareMesh &Objs type = (
	
	
		for o = 1 to Objs.count do
			(	
					if PBRT_Render.exportObjs.checked then (
						
							(
								
									MeshMod = Turn_To_Mesh()()				
									MeshMod.name = "Mesh_Edit"
								
									if type == "instance" then  (
										
								--		if classof Objs[o][1] != Sphere then addModifier Objs[o][1] MeshMod
										addModifier Objs[o][1] MeshMod
									) else	(
	--									if classof Objs[o] != Sphere then addModifier Objs[o] MeshMod	
										addModifier Objs[o] MeshMod
									)
								
						)
												
					)							
				
			)
	
	
)

-- function to remove the mesh modofier
fn cleanMesh &Objs type = (
	
		for o = 1 to Objs.count do
			(				
									
				if PBRT_Render.exportObjs.checked then (
						if type == "instance" then (
--									if classof Objs[o][1] != Sphere and then deleteModifier  Objs[o][1] Objs[o][1].Mesh_Edit
											deleteModifier  Objs[o][1] Objs[o][1].Mesh_Edit
									)else (
	--										if classof Objs[o] != Sphere then deleteModifier  Objs[o][1] Objs[o][1].Mesh_Edit
											deleteModifier  Objs[o][1] Objs[o][1].Mesh_Edit
											)
							)	
											
				
			)
		
)

-- function that return an array of subarrays holding the instances
-- eje arr[a][b] a is the main object and b the instances of a
fn getInstanceArray TheObjects = (
	
				Local InstanceObjs = #()
				Local doneObjs = #() 
				
				-- sort the objects Instances

			for o = 1 to TheObjects.count where (findItem doneObjs TheObjects[o] == 0 ) do (
					instanceArr = #()

		--			format "Master Obj:% \n" TheObjects[o].name
					
				   instObjs = #()
				   InstanceMgr.GetInstances TheObjects[o] &instObjs -- find the instances of the first obj

				   for thisInst in instObjs do (

		--			format "MasObj:% Instance:% \n" TheObjects[o].name thisInst.name 
					 append doneObjs thisInst -- mark this object done	
					   
					if (not thisInst.isHidden) do append instanceArr thisInst -- append the instance to the instanarra
				   )

				   append InstanceObjs instanceArr -- append the group of instances to the main array
				)
				
				
		return InstanceObjs
)

fn Pbrt_Export = (
			
	Local dirPath = PBRT_System.PBRT_SavePath_edt.text 
	Only = filterstring dirPath "\\"
	newTxt = substring dirPath 1 (dirPath.count - only[only.count].count - 1)
	dirPath = newTxt + "\\"
	Local theFileOnly = filterstring PBRT_System.PBRT_SavePath_edt.text "\\"		
	
	-- we need to check if the user provided file extension so we dont duplicat it
	Local theName = theFileOnly[theFileOnly.count] 	-- name pat of the complete string with or without extension, previously we already spli it by \
	Local remExt =  filterstring theName "."	-- split it again now by .

	-- if it has cont more then 1 it means it has an extension
	if remExt.count > 1 then theFileOnly = remExt[1] -- assign the element 1, removing the extension
		else theFileOnly = remExt[remExt.count] 		-- the user did not provide the extension
		
	-- now we make the complete names.
	Local mainTexFile = dirPath + theFileOnly + ".pbrt"
	Local geomTexFile = dirPath + theFileOnly+ "-geom" + ".pbrt"
	Local matTexFile = dirPath + theFileOnly + "-mat" + ".pbrt"

	Local out_Main_file = createfile mainTexFile
	
	if PBRT_Render.exportObjs.checked then Local out_Geom_file = createfile geomTexFile
	if PBRT_Render.exportMaps.checked then Local out_Mat_file = createfile  matTexFile
	
-----------------------------------------------------------------------------------------------------------------------------------------------------	

	format "# Gtool Exporter Ver. % for pbrt v2.0 Render, www.pbrt.org \n" ver to:out_Main_file
	format "# Created by Guillermo M. Leal LLaguno,  g.llaguno@evvisual.com,  www.evvisual.com \n"  to:out_Main_file
	format "# Exported on: % \n\n"  localTime to:out_Main_file
	
	format "Scale -1 1 1  # account for fixed lookat bug... \n\n" to:out_Main_file
	
	Pbrt_ExpCamera CurrentCamera	out_Main_file	
	Pbrt_PixelFilter PBRT_Render.Render_Filter_ddl.selected	out_Main_file -- FilterType[3] -- select the type of filter from the filtertype array
	Pbrt_Sampler PBRT_Render.Render_Sampler_ddl.selection	out_Main_file   --SamplerType[1] -- select the type of sampler from the SamplerType array
	Pbrt_SurfaceIntegrator PBRT_Render.Render_Integ_ddl.selection	out_Main_file  --integratorType[1] -- select the type of sampler from the SamplerType array
	Pbrt_VolumeIntegrator PBRT_Render.Render_Vol_ddl.selection	out_Main_file --VolumeType[1] -- select the type of sampler from the SamplerType array
	Pbrt_Renderer PBRT_Render.Render_Rnd_ddl.selection	out_Main_file --RenderTypes[2] -- select the type of sampler from the SamplerType array
	Pbrt_Accelerator PBRT_System.Accelerator_ddl.selection out_Main_file


	format "\n\n" 	to:out_Main_file
	format "WorldBegin\n" 	to:out_Main_file
	format "\n" 	to:out_Main_file


	Pbrt_Enviorment	out_Main_file

	-- Check for Groups and if it finds them, open them
	local TheGroups = #()
		
	for z in objects do 
		(
		if (isGroupMember z) then
			(
				setGroupMember z false
				append TheGroups z
			)
		)
	
	TheLights 	=  for o in lights where  classof o != Targetobject  collect o
	TheObjects = for o in objects where superclassof o == GeometryClass and classof o != Targetobject and (not o.isHidden) collect o
	
	-- Lights
	format "# ***** Lights ***** \n\n" to:out_Main_file
	for o = 1 to TheLights.count do
		(
			-- check if the light is on
			if TheLights[o].enabled == on then (
				Pbrt_Lights TheLights[o] out_Main_file
				)
		)
	format "# ***** End  Lights *****\n" to:out_Main_file
	
	-- ***** NOTE ****
/*
the instances have been a paint to do!!!!!
for now they are working but i dont know why i have to recreatethe instance array and add a convert to mesh modifier??????
*/		
		Local InstanceObjs = getInstanceArray TheObjects

		if PBRT_System.PBRT_UseInstance_chk.checked and PBRT_Render.exportObjs.checked do (prepareMesh InstanceObjs  "instance")
				
		-- check if we want to use instances 
		if PBRT_System.PBRT_UseInstance_chk.checked then (

				Local InstanceObjs2 = #()
				Local doneObjs = #() 
				
				-- sort the objects Instances

			for o = 1 to TheObjects.count where (findItem doneObjs TheObjects[o] == 0 ) do (
					instanceArr = #()

		--			format "Master Obj:% \n" TheObjects[o].name
					
				   instObjs = #()
				   InstanceMgr.GetInstances TheObjects[o] &instObjs -- find the instances of the first obj

				   for thisInst in instObjs do (

		--			format "MasObj:% Instance:% \n" TheObjects[o].name thisInst.name 
					 append doneObjs thisInst -- mark this object done	
					if (not thisInst.isHidden) do append instanceArr thisInst -- append the instance to the instanarra
				   )

				   append InstanceObjs2 instanceArr -- append the group of instances to the main array
				)
				
					
			-- now procees the objects
			for o = 1 to InstanceObjs2.count do
					(						
--						format "%\n" InstanceObjs[o].count
						-- if there is more then one in the array then they are instances
						if InstanceObjs2[o].count > 1 then (
							if PBRT_Render.exportMaps.checked then Pbrt_Material InstanceObjs2[o][1] out_Mat_file	
										
							if PBRT_Render.exportObjs.checked then Pbrt_GeometryInstance InstanceObjs2[o] out_Geom_file
							
						)else -- Single Objects
								(
								--	format "%\n" (classof InstanceObjs[o][1])
										 -- check if we want to export the maps, maybe they are already exported and we are just chageing render setings
									if PBRT_Render.exportMaps.checked then (
										 Pbrt_Material InstanceObjs2[o][1] out_Mat_file	
										)	
										
									-- check if we want to exprt geometry, in case we are re rendering with new settings
									if PBRT_Render.exportObjs.checked then (
										Pbrt_Geometry InstanceObjs2[o][1] out_Geom_file
										)								
								)
					-- nesecito modicifar esto paa que incluya el proceso de cada instance
					PBRT_Render.export_prog.value = 100.000 / InstanceObjs2.count * o
					progress_text = ((100.0 / InstanceObjs2.count * o) as integer)as string  + " % done"
					pushprompt progress_text	
								
					)
				PBRT_Render.export_prog.value = 0
				)
				else (

					prepareMesh TheObjects "single"
				
						for o = 1 to TheObjects.count do
								(			
										 -- check if we want to export the maps, maybe they are already exported and we are just chageing render setings
										if PBRT_Render.exportMaps.checked then (
											 Pbrt_Material TheObjects[o] out_Mat_file	
											)	
											
										-- check if we want to exprt geometry, in case we are re rendering with new settings
										if PBRT_Render.exportObjs.checked then (
											Pbrt_Geometry TheObjects[o] out_Geom_file
											)
									
									PBRT_Render.export_prog.value = 100.000 / TheObjects.count * o
									progress_text = ((100.0 / TheObjects.count * o) as integer)as string  + " % done"
									pushprompt progress_text									
								)
							PBRT_Render.export_prog.value = 0
					cleanMesh  InstanceObjs "single"		
				)
		
		
		if PBRT_System.PBRT_UseInstance_chk.checked and PBRT_Render.exportObjs.checked do (cleanMesh  InstanceObjs "instance")

	-- if it was in a group put it back
	for z in TheGroups do setGroupMember z true
	
	format "Include \"%\" \n" (theFileOnly+ "-mat" + ".pbrt") 	to:out_Main_file
	format "Include \"%\" \n" (theFileOnly + "-geom" + ".pbrt") 	to:out_Main_file

	format "\nWorldEnd\n" 	to:out_Main_file
		
	------------------------------------------------------------------------------------------------
		
	close out_Main_file 

	if PBRT_Render.exportObjs.checked then close out_Geom_file
	if PBRT_Render.exportMaps.checked then close out_Mat_file
	
	)
	
--- **************************************************************************************************************************************
--- ****************  GUI *********************
	
rollout PBRT_System "System" width:304 height:233
	(
		label 				PBRT_Save_lbl 	"Save path:" 			pos:[5,9] width:56 height:17
		editText 			PBRT_SavePath_edt 	"" 							pos:[57,8] width:165 height:17
		dropdownList 	PBRT_Ext_ddl 		"" 							pos:[225,7] width:48 height:17 items:#(".exr", ".tga")	selection:1	
		button 				PBRT_Save_btn 	"..." 						pos:[278,8] width:15 height:17
		
		label 				PBRT_exe_lbl 		"PBRT exe:" 			pos:[4,31] width:94 height:17
		editText 			PBRT_exe_Path_edt		"" 							pos:[97,30] width:174 height:17	
		button 				PBRT_exe_Path_btn 	"..." 						pos:[278,30] width:15 height:17
		
		label 				Threads_lbl 		"Threads:" 				pos:[4,54] width:47 height:17	
		checkbutton		ThreadsAuto_btn "Auto Detect" 			pos:[52,53] width:64 height:17	checked:true
		spinner 			Threads_spn 		"" 							pos:[124,54] width:65 height:16 range:[0,100,1] type:#integer

		checkbox 			PBRT_UseInstance_chk "Use Instances" pos:[200,54] width:100 height:17 checked:true
		spinner 			PBRT_WriteFrequency 		"Write Frequency" pos:[40,74] width:120 height:16 range:[-1,1000000,-1] type:#integer
		checkbox 			PBRT_DisableMaps_chk "Disable Maps" pos:[200,74] width:100 height:17 
		
		groupBox 			Accelerator_grp "Accelerator" pos:[3,90] width:292 height:115
		dropDownList 	Accelerator_ddl "" pos:[11,106] width:66 height:21 items:#("None", "Grid", "KdTree", "Bvh") selection:4

		checkbutton 		Acce_Grid_Refine_btn "Refine Immediately" pos:[99,108] width:100 height:17	 	visible:false
		
		label  			Acce_KdTree_Intersect_lbl 		"Intersect Cost:"  pos:[95,105] 	visible:false width:78 height:16 
		spinner 		Acce_KdTree_Intersect_spn 		"" 						pos:[172,104] visible:false width:60 height:16 range:[0,1000,80] type:#integer
		label 			Acce_KdTree_Traversal_lbl		"Traversal Cost:" 	pos:[95,125] 	visible:false width:78 height:16
		spinner 		Acce_KdTree_Traversal_spn 		"" 						pos:[172,124] visible:false width:60 height:16 range:[0,1000,1] type:#integer
		label 			Acce_KdTree_EmptyBonus_lbl 	"Empty Bonus:" 	pos:[96,144] 	visible:false width:78 height:16
		spinner 		Acce_KdTree_EmptyBonus_Spn 	"" 						pos:[172,143] visible:false width:60 height:16 range:[0,100,0.2] type:#float
		label 			Acce_KdTree_MaxPrims_lbl 		"Max Primes:" 		pos:[95,163] 	visible:false width:78 height:16
		spinner 		Acce_KdTree_MaxPrimes_spn 	"" 						pos:[172,162] visible:false width:60 height:16 range:[0,1000,1] type:#integer
		label 			Acce_KdTree_MaxDepth_lbl 		"Max Depth:" 		pos:[96,182] 	visible:false width:78 height:16
		spinner 		Acce_KdTree_MaxDepth_spn 		"" 						pos:[173,181] visible:false width:60 height:16 range:[0,100,-1] type:#integer

		button 			SettingSave_btn_ "Save Settings" pos:[5,210] width:96 height:17
		button 			SettingsLoad_btn "Load Settings" pos:[106,210] width:96 height:17
		
		checkbox 		usePrimitives_chk "Use Primitives" pos:[5,230] width:100 height:17	 	visible:true

		
		on PBRT_System open do (
			PBRT_exe_Path_edt.text 		= PBRT_EXE_PATH
			PBRT_SavePath_edt.text 		= SAVEPATH
		)
		
		on PBRT_Save_btn pressed do
			(		
				SAVEPATH = getSaveFileName caption:"-- Gtools --  Select File Name and Render Path. NO extension" types:"PBRT(*.pbrt)|*.pbrt|All|*.*|"
				if (SAVEPATH  != undefined) then
				(			
				--	OrifilePathReduced =  substring SAVEPATH 1 -1 -- reduces the path to fit on the window
					PBRT_SavePath_edt.text = SAVEPATH --"\.. " + (OrifilePathReduced as string)	
				)		
			)
			
		on PBRT_exe_Path_btn pressed do
			(		
				PBRT_EXE_PATH = getSavePath caption:"-- Gtools --  Select the Location of Pbrt.exe"
				if (PBRT_EXE_PATH  != undefined) then
				(			
					--OrifilePathReduced =  substring PBRT_EXE_PATH 1 -1 -- reduces the path to fit on the window
					PBRT_exe_Path_edt.text = PBRT_EXE_PATH-- "\.. " + (OrifilePathReduced as string)	
				)		
			)
			
		on Accelerator_ddl selected val do (
			case Accelerator_ddl.selection of 
				(
					1:
					(
					Acce_Grid_Refine_btn.visible = false
					Acce_KdTree_Intersect_lbl.visible = false
					Acce_KdTree_Intersect_spn.visible = false
					Acce_KdTree_Traversal_lbl.visible = false
					Acce_KdTree_Traversal_spn.visible = false
					Acce_KdTree_EmptyBonus_lbl.visible = false
					Acce_KdTree_EmptyBonus_Spn.visible = false
					Acce_KdTree_MaxPrims_lbl.visible = false
					Acce_KdTree_MaxPrimes_spn.visible = false
					Acce_KdTree_MaxDepth_lbl.visible = false
					Acce_KdTree_MaxDepth_spn.visible = false
					)
					2:
					(
					Acce_Grid_Refine_btn.visible = true					
					Acce_KdTree_Intersect_lbl.visible = false
					Acce_KdTree_Intersect_spn.visible = false
					Acce_KdTree_Traversal_lbl.visible = false
					Acce_KdTree_Traversal_spn.visible = false
					Acce_KdTree_EmptyBonus_lbl.visible = false
					Acce_KdTree_EmptyBonus_Spn.visible = false
					Acce_KdTree_MaxPrims_lbl.visible = false
					Acce_KdTree_MaxPrimes_spn.visible = false
					Acce_KdTree_MaxDepth_lbl.visible = false
					Acce_KdTree_MaxDepth_spn.visible = false
					)				
					3:
					(
					Acce_Grid_Refine_btn.visible = false					
					Acce_KdTree_Intersect_lbl.visible = true
					Acce_KdTree_Intersect_spn.visible = true
					Acce_KdTree_Traversal_lbl.visible = true
					Acce_KdTree_Traversal_spn.visible = true
					Acce_KdTree_EmptyBonus_lbl.visible = true
					Acce_KdTree_EmptyBonus_Spn.visible = true
					Acce_KdTree_MaxPrims_lbl.visible = true
					Acce_KdTree_MaxPrimes_spn.visible = true
					Acce_KdTree_MaxDepth_lbl.visible = true
					Acce_KdTree_MaxDepth_spn.visible = true
					)				
					4:
					(
					Acce_Grid_Refine_btn.visible = false					
					Acce_KdTree_Intersect_lbl.visible = false
					Acce_KdTree_Intersect_spn.visible = false
					Acce_KdTree_Traversal_lbl.visible = false
					Acce_KdTree_Traversal_spn.visible = false
					Acce_KdTree_EmptyBonus_lbl.visible = false
					Acce_KdTree_EmptyBonus_Spn.visible = false
					Acce_KdTree_MaxPrims_lbl.visible = false
					Acce_KdTree_MaxPrimes_spn.visible = false
					Acce_KdTree_MaxDepth_lbl.visible = false
					Acce_KdTree_MaxDepth_spn.visible = false
					)				
				)
		)
		
	)

rollout PBRT_Output "Output" 
	(

		
		fn  load_icon  = 
			 (
				a = openbitmap  ((GetDir #ui)+ "\\Icons\\"  + "SnapsControls_16i.bmp")
				return a
			 )
		 
		spinner RenderWidth_spn "Render Width:" 	pos:[55,10] width:105 height:16 range:[0,100000,0] type:#integer
		spinner RenderHeight_spn "Render Height:" pos:[55,30] width:105 height:16 range:[0,100000,0] type:#integer
		checkbutton RenderLock "Button" pos:[170,18] width:17 height:17 images:#(load_icon (), undefined, 38, 1, 3, 3, 3)
		checkbutton crop_btn "Use Crop" pos:[210,19] width:70 height:16
		
		on PBRT_Output open do (
			RenderWidth_spn.value = renderwidth
			RenderHeight_spn.value = renderheight
			ResRatio = (RenderWidth_spn.value as float) / (RenderHeight_spn.value  as float)
		)
		
		on RenderWidth_spn changed val do (
			if RenderLock.checked do RenderHeight_spn.value = ((val / ResRatio) as integer)
			renderWidth = val
			renderHeight = RenderHeight_spn.value	
		)
		
		on RenderHeight_spn changed val do(
		--	ResRatio = (RenderWidth_spn.value as float) / (RenderHeight_spn.value  as float)
			
			if RenderLock.checked do RenderWidth_spn.value = ((val / ResRatio) as integer)
			renderWidth = RenderWidth_spn.value
			renderHeight = val
		)
		
	)
 
fn PBRT_SaveProperties = (
	
	-- if there is no helper create it
	if $PBRT_Helper == undefined then (
			thePoint = Point pos:[0,0,0] isSelected:off centermarker:on cross:off axistripod:off Box:off
			thePoint.name = "PBRT_Helper"
		)
		
	-- Set the Properties --
		
	-- **** RENDER SETTINGS **** --
		-- Sampler
	setuserprop $PBRT_Helper "Render_Sampler_ddl" PBRT_Render.Render_Sampler_ddl.selection
	setuserprop $PBRT_Helper "Render_Sampler_spn" PBRT_Render.Render_Sampler_spn.value
		-- Render
	setuserprop $PBRT_Helper "Render_Rnd_ddl" PBRT_Render.Render_Rnd_ddl.selection		
	setuserprop $PBRT_Helper "Render_Rnd_Aggre_Niters_spn" PBRT_Render.Render_Rnd_Aggre_Niters_spn.value		

	setuserprop $PBRT_Helper "Render_Rnd_Probe_Lmax_spn" PBRT_Render.Render_Rnd_Probe_Lmax_spn.value	
	setuserprop $PBRT_Helper "Render_Rnd_Probe_MaxDepth_spn" PBRT_Render.Render_Rnd_Probe_MaxDepth_spn.value
	setuserprop $PBRT_Helper "Render_Rnd_Probe_IndSamp_spn" PBRT_Render.Render_Rnd_Probe_IndSamp_spn.value
	setuserprop $PBRT_Helper "Render_Rnd_Probe_SampSpace_spn" PBRT_Render.Render_Rnd_Probe_SampSpace_spn.value	
	setuserprop $PBRT_Helper "Render_Rnd_Probe_Direct_btn" PBRT_Render.Render_Rnd_Probe_Direct_btn.state
	setuserprop $PBRT_Helper "Render_Rnd_Probe_Indirect_btn" PBRT_Render.Render_Rnd_Probe_Indirect_btn.state	
	setuserprop $PBRT_Helper "Render_Rnd_Probe_time_spn" PBRT_Render.Render_Rnd_Probe_time_spn.value		
	
	setuserprop $PBRT_Helper "Render_Rnd_Metro_Probab_spn" PBRT_Render.Render_Rnd_Metro_Probab_spn.value
	setuserprop $PBRT_Helper "Render_Rnd_Metro_Nsamples_spn" PBRT_Render.Render_Rnd_Metro_Nsamples_spn.value
	setuserprop $PBRT_Helper "Render_Rnd_Metro_Samples_spn" PBRT_Render.Render_Rnd_Metro_Samples_spn.value
	setuserprop $PBRT_Helper "Render_Rnd_Metro_BootStr_spn" PBRT_Render.Render_Rnd_Metro_BootStr_spn.value
	setuserprop $PBRT_Helper "Render_Rnd_Metro_DirectSamp" PBRT_Render.Render_Rnd_Metro_DirectSamp.value		
	setuserprop $PBRT_Helper "Render_Rnd_Metro_DoDirect" PBRT_Render.Render_Rnd_Metro_DoDirect.state	
	setuserprop $PBRT_Helper "Render_Rnd_Metro_MaxReje_spn" PBRT_Render.Render_Rnd_Metro_MaxReje_spn.value
	setuserprop $PBRT_Helper "Render_Rnd_Metro_MaxDepth" PBRT_Render.Render_Rnd_Metro_MaxDepth.value
	setuserprop $PBRT_Helper "Render_Rnd_Metro_Indirect_btn" PBRT_Render.Render_Rnd_Metro_Indirect_btn.state
	
	setuserprop $PBRT_Helper "Render_Rnd_Surface_MinSam_spn" PBRT_Render.Render_Rnd_Surface_MinSam_spn.value
		
	setuserprop $PBRT_Helper "Render_Integ_ddl" PBRT_Render.Render_Integ_ddl.selection
	
	setuserprop $PBRT_Helper "Render_Integ_AO_Samples_spn" PBRT_Render.Render_Integ_AO_Samples_spn.value
	setuserprop $PBRT_Helper "Render_Integ_AO_MaxDist_spn" PBRT_Render.Render_Integ_AO_MaxDist_spn.value
	
	setuserprop $PBRT_Helper "Render_Integ_Diffuse_Max_spn" PBRT_Render.Render_Integ_Diffuse_Max_spn.value
	setuserprop $PBRT_Helper "Render_Integ_Diffuse_Ns_spn" PBRT_Render.Render_Integ_Diffuse_Ns_spn.value
	
	setuserprop $PBRT_Helper "Render_Integ_Dipole_MaxDepth_spn" PBRT_Render.Render_Integ_Dipole_MaxDepth_spn.value
	setuserprop $PBRT_Helper "Render_Integ_Dipole_MaxError_spn" PBRT_Render.Render_Integ_Dipole_MaxError_spn.value
	setuserprop $PBRT_Helper "Render_Integ_Dipole_MinDist_spn" PBRT_Render.Render_Integ_Dipole_MinDist_spn.value
	
	setuserprop $PBRT_Helper "Render_Integ_Direct_Depth_spn" PBRT_Render.Render_Integ_Direct_Depth_spn.value
	
	setuserprop $PBRT_Helper "Render_Integ_Emission_Step_spn" PBRT_Render.Render_Integ_Emission_Step_spn.value
	
	setuserprop $PBRT_Helper "Render_Integ_Glossy_LMax_spn" PBRT_Render.Render_Integ_Glossy_LMax_spn.value
	setuserprop $PBRT_Helper "Render_Integ_Glossy_Ns_spn" PBRT_Render.Render_Integ_Glossy_Ns_spn.value
	setuserprop $PBRT_Helper "Render_Integ_Glossy_DotR_btn" PBRT_Render.Render_Integ_Glossy_DotR_btn.state
	setuserprop $PBRT_Helper "Render_Integ_Glossy_Kd_spn" PBRT_Render.Render_Integ_Glossy_Kd_spn.value
	setuserprop $PBRT_Helper "Render_Integ_Glossy_Ks_spn" PBRT_Render.Render_Integ_Glossy_Ks_spn.value
	setuserprop $PBRT_Helper "Render_Integ_Glossy_Rough_spn" PBRT_Render.Render_Integ_Glossy_Rough_spn.value
	
	setuserprop $PBRT_Helper "Render_Integ_IGI_LightPath_spn" PBRT_Render.Render_Integ_IGI_LightPath_spn.value
	setuserprop $PBRT_Helper "Render_Integ_IGI_LightSet_spn" PBRT_Render.Render_Integ_IGI_LightSet_spn.value
	setuserprop $PBRT_Helper "Render_Integ_IGI_MinDist_spn" PBRT_Render.Render_Integ_IGI_MinDist_spn.value
	setuserprop $PBRT_Helper "Render_Integ_IGI_RRThrs_spn" PBRT_Render.Render_Integ_IGI_RRThrs_spn.value
	setuserprop $PBRT_Helper "Render_Integ_IGI_MaxDepth_spn" PBRT_Render.Render_Integ_IGI_MaxDepth_spn.value
	setuserprop $PBRT_Helper "Render_Integ_IGI_GLimit_spn" PBRT_Render.Render_Integ_IGI_GLimit_spn.value
	setuserprop $PBRT_Helper "Render_Integ_IGI_Gather_spn" PBRT_Render.Render_Integ_IGI_Gather_spn.value
	
	setuserprop $PBRT_Helper "Render_Integ_IrCach_MinWht_spn" PBRT_Render.Render_Integ_IrCach_MinWht_spn.value
	setuserprop $PBRT_Helper "Render_Integ_IrCach_MinSpac_spn" PBRT_Render.Render_Integ_IrCach_MinSpac_spn.value
	setuserprop $PBRT_Helper "Render_Integ_IrCach_MaxSpac_spn" PBRT_Render.Render_Integ_IrCach_MaxSpac_spn.value
	setuserprop $PBRT_Helper "Render_Integ_IrCach_MaxAngle_spn" PBRT_Render.Render_Integ_IrCach_MaxAngle_spn.value
	setuserprop $PBRT_Helper "Render_Integ_IrCach_SpecDepth_spn" PBRT_Render.Render_Integ_IrCach_SpecDepth_spn.value
	setuserprop $PBRT_Helper "Render_Integ_IrCach_indDepth_spn" PBRT_Render.Render_Integ_IrCach_indDepth_spn.value
	setuserprop $PBRT_Helper "Render_Integ_IrCach_Samples_spn" PBRT_Render.Render_Integ_IrCach_Samples_spn.value
	
	setuserprop $PBRT_Helper "Render_Integ_Path_Depth_spn" PBRT_Render.Render_Integ_Path_Depth_spn.value
	
	setuserprop $PBRT_Helper "Render_Integ_Photon_Caustic_spn" PBRT_Render.Render_Integ_Photon_Caustic_spn.value
	setuserprop $PBRT_Helper "Render_Integ_Photon_Indirect_spn" PBRT_Render.Render_Integ_Photon_Indirect_spn.value
	setuserprop $PBRT_Helper "Render_Integ_Photon_nUsed_spn" PBRT_Render.Render_Integ_Photon_nUsed_spn.value
	setuserprop $PBRT_Helper "Render_Integ_Photon_SpecDepth_spn" PBRT_Render.Render_Integ_Photon_SpecDepth_spn.value
	setuserprop $PBRT_Helper "Render_Integ_Photon_MaxPhoDepth_spn" PBRT_Render.Render_Integ_Photon_MaxPhoDepth_spn.value
	setuserprop $PBRT_Helper "Render_Integ_Photon_FinalGath_btn" PBRT_Render.Render_Integ_Photon_FinalGath_btn.state
	setuserprop $PBRT_Helper "Render_Integ_Photon_gatherSamp_spn" PBRT_Render.Render_Integ_Photon_gatherSamp_spn.value
	setuserprop $PBRT_Helper "Render_Integ_Photon_MaxDist_spn" PBRT_Render.Render_Integ_Photon_MaxDist_spn.value
	setuserprop $PBRT_Helper "Render_Integ_Photon_GathAng_spn" PBRT_Render.Render_Integ_Photon_GathAng_spn.value
	
	setuserprop $PBRT_Helper "Render_Integ_Single_Step_spn" PBRT_Render.Render_Integ_Single_Step_spn.value
	
	setuserprop $PBRT_Helper "Render_Integ_Whitted_MaDepth_spn" PBRT_Render.Render_Integ_Whitted_MaDepth_spn.value
	
	setuserprop $PBRT_Helper "Render_Vol_ddl" PBRT_Render.Render_Vol_ddl.selection
	setuserprop $PBRT_Helper "Render_Vol_Step_spn" PBRT_Render.Render_Vol_Step_spn.value
	
	setuserprop $PBRT_Helper "Render_Filter_ddl" PBRT_Render.Render_Filter_ddl.selection
	
	setuserprop $PBRT_Helper "Render_Filter_Box_X" PBRT_Render.Render_Filter_Box_X.value
	setuserprop $PBRT_Helper "Render_Filter_Box_Y" PBRT_Render.Render_Filter_Box_Y.value
	
	setuserprop $PBRT_Helper "Render_Filter_Gauss_X" PBRT_Render.Render_Filter_Gauss_X.value
	setuserprop $PBRT_Helper "Render_Filter_Gauss_Y" PBRT_Render.Render_Filter_Gauss_Y.value
	setuserprop $PBRT_Helper "Render_Filter_Gauss_Alpah" PBRT_Render.Render_Filter_Gauss_Alpah.value
	
	setuserprop $PBRT_Helper "Render_Filter_Mitch_X" PBRT_Render.Render_Filter_Mitch_X.value
	setuserprop $PBRT_Helper "Render_Filter_Mitch_Y" PBRT_Render.Render_Filter_Mitch_Y.value
	setuserprop $PBRT_Helper "Render_Filter_Mitch_Sharp" PBRT_Render.Render_Filter_Mitch_Sharp.value
	
	setuserprop $PBRT_Helper "Render_Filter_Sinc_X" PBRT_Render.Render_Filter_Sinc_X.value
	setuserprop $PBRT_Helper "Render_Filter_Sinc_Y" PBRT_Render.Render_Filter_Sinc_Y.value
	setuserprop $PBRT_Helper "Render_Filter_Sinc_Tau" PBRT_Render.Render_Filter_Sinc_Tau.value
	
	setuserprop $PBRT_Helper "Render_Filter_Tri_X" PBRT_Render.Render_Filter_Tri_X.value
	setuserprop $PBRT_Helper "Render_Filter_Tri_Y" PBRT_Render.Render_Filter_Tri_Y.value

	setuserprop $PBRT_Helper "exportObjs" PBRT_Render.exportObjs.state  
	setuserprop $PBRT_Helper "exportMaps" PBRT_Render.exportMaps.state 
	setuserprop $PBRT_Helper "envMult" PBRT_Render.envMult.value	
	setuserprop $PBRT_Helper "envSamp" PBRT_Render.envSamp.value
	
	setuserprop $PBRT_Helper "useGlobalMat_chk" PBRT_Render.useGlobalMat_chk.state
		
		
)

fn PBRT_GetProperties = (
	
	-- if there is no helper create it
	if $PBRT_Helper != undefined then (
			-- Get the Properties --

			-- **** RENDER SETTINGS **** --
				-- Sampler
			PBRT_Render.Render_Sampler_ddl.selection = getuserprop $PBRT_Helper "Render_Sampler_ddl" 
			PBRT_Render.Render_Sampler_spn.value = getuserprop $PBRT_Helper "Render_Sampler_spn" 
				-- Render
			PBRT_Render.Render_Rnd_ddl.selection = getuserprop $PBRT_Helper "Render_Rnd_ddl" 		
			PBRT_Render.Render_Rnd_Aggre_Niters_spn.value = getuserprop $PBRT_Helper "Render_Rnd_Aggre_Niters_spn" 		

			PBRT_Render.Render_Rnd_Probe_Lmax_spn.value =getuserprop $PBRT_Helper "Render_Rnd_Probe_Lmax_spn"  	
			PBRT_Render.Render_Rnd_Probe_MaxDepth_spn.value = getuserprop $PBRT_Helper "Render_Rnd_Probe_MaxDepth_spn" 
			PBRT_Render.Render_Rnd_Probe_IndSamp_spn.value = getuserprop $PBRT_Helper "Render_Rnd_Probe_IndSamp_spn" 
			PBRT_Render.Render_Rnd_Probe_SampSpace_spn.value = getuserprop $PBRT_Helper "Render_Rnd_Probe_SampSpace_spn"	
			PBRT_Render.Render_Rnd_Probe_Direct_btn.state = getuserprop $PBRT_Helper "Render_Rnd_Probe_Direct_btn" 
			PBRT_Render.Render_Rnd_Probe_Indirect_btn.state = getuserprop $PBRT_Helper "Render_Rnd_Probe_Indirect_btn" 
			PBRT_Render.Render_Rnd_Probe_time_spn.value = getuserprop $PBRT_Helper "Render_Rnd_Probe_time_spn" 		
			
			PBRT_Render.Render_Rnd_Metro_Probab_spn.value = getuserprop $PBRT_Helper "Render_Rnd_Metro_Probab_spn" 
			PBRT_Render.Render_Rnd_Metro_Nsamples_spn.value = getuserprop $PBRT_Helper "Render_Rnd_Metro_Nsamples_spn" 
			PBRT_Render.Render_Rnd_Metro_Samples_spn.value = getuserprop $PBRT_Helper "Render_Rnd_Metro_Samples_spn" 
			PBRT_Render.Render_Rnd_Metro_BootStr_spn.value = getuserprop $PBRT_Helper "Render_Rnd_Metro_BootStr_spn" 
			PBRT_Render.Render_Rnd_Metro_DirectSamp.value = getuserprop $PBRT_Helper "Render_Rnd_Metro_DirectSamp"		
			PBRT_Render.Render_Rnd_Metro_DoDirect.state = getuserprop $PBRT_Helper "Render_Rnd_Metro_DoDirect" 	
			PBRT_Render.Render_Rnd_Metro_MaxReje_spn.value = getuserprop $PBRT_Helper "Render_Rnd_Metro_MaxReje_spn" 
			PBRT_Render.Render_Rnd_Metro_MaxDepth.value = getuserprop $PBRT_Helper "Render_Rnd_Metro_MaxDepth" 
			PBRT_Render.Render_Rnd_Metro_Indirect_btn.state = getuserprop $PBRT_Helper "Render_Rnd_Metro_Indirect_btn" 
			
			PBRT_Render.Render_Rnd_Surface_MinSam_spn.value = getuserprop $PBRT_Helper "Render_Rnd_Surface_MinSam_spn" 
				
			PBRT_Render.Render_Integ_ddl.selection = getuserprop $PBRT_Helper "Render_Integ_ddl" 
			
			PBRT_Render.Render_Integ_AO_Samples_spn.value = getuserprop $PBRT_Helper "Render_Integ_AO_Samples_spn" 
			PBRT_Render.Render_Integ_AO_MaxDist_spn.value = getuserprop $PBRT_Helper "Render_Integ_AO_MaxDist_spn" 
			
			PBRT_Render.Render_Integ_Diffuse_Max_spn.value = getuserprop $PBRT_Helper "Render_Integ_Diffuse_Max_spn" 
			PBRT_Render.Render_Integ_Diffuse_Ns_spn.value = getuserprop $PBRT_Helper "Render_Integ_Diffuse_Ns_spn"
			
			PBRT_Render.Render_Integ_Dipole_MaxDepth_spn.value = getuserprop $PBRT_Helper "Render_Integ_Dipole_MaxDepth_spn" 
			PBRT_Render.Render_Integ_Dipole_MaxError_spn.value = getuserprop $PBRT_Helper "Render_Integ_Dipole_MaxError_spn" 
			PBRT_Render.Render_Integ_Dipole_MinDist_spn.value = getuserprop $PBRT_Helper "Render_Integ_Dipole_MinDist_spn" 
			
			PBRT_Render.Render_Integ_Direct_Depth_spn.value = getuserprop $PBRT_Helper "Render_Integ_Direct_Depth_spn" 
			
			PBRT_Render.Render_Integ_Emission_Step_spn.value = getuserprop $PBRT_Helper "Render_Integ_Emission_Step_spn" 
			
			PBRT_Render.Render_Integ_Glossy_LMax_spn.value = getuserprop $PBRT_Helper "Render_Integ_Glossy_LMax_spn" 
			PBRT_Render.Render_Integ_Glossy_Ns_spn.value = getuserprop $PBRT_Helper "Render_Integ_Glossy_Ns_spn" 
			PBRT_Render.Render_Integ_Glossy_DotR_btn.state = getuserprop $PBRT_Helper "Render_Integ_Glossy_DotR_btn" 
			PBRT_Render.Render_Integ_Glossy_Kd_spn.value = getuserprop $PBRT_Helper "Render_Integ_Glossy_Kd_spn" 
			PBRT_Render.Render_Integ_Glossy_Ks_spn.value = getuserprop $PBRT_Helper "Render_Integ_Glossy_Ks_spn" 
			PBRT_Render.Render_Integ_Glossy_Rough_spn.value = getuserprop $PBRT_Helper "Render_Integ_Glossy_Rough_spn" 
			
			PBRT_Render.Render_Integ_IGI_LightPath_spn.value = getuserprop $PBRT_Helper "Render_Integ_IGI_LightPath_spn" 
			PBRT_Render.Render_Integ_IGI_LightSet_spn.value = getuserprop $PBRT_Helper "Render_Integ_IGI_LightSet_spn" 
			PBRT_Render.Render_Integ_IGI_MinDist_spn.value = getuserprop $PBRT_Helper "Render_Integ_IGI_MinDist_spn" 
			PBRT_Render.Render_Integ_IGI_RRThrs_spn.value = getuserprop $PBRT_Helper "Render_Integ_IGI_RRThrs_spn" 
			PBRT_Render.Render_Integ_IGI_MaxDepth_spn.value = getuserprop $PBRT_Helper "Render_Integ_IGI_MaxDepth_spn" 
			PBRT_Render.Render_Integ_IGI_GLimit_spn.value = getuserprop $PBRT_Helper "Render_Integ_IGI_GLimit_spn" 
			PBRT_Render.Render_Integ_IGI_Gather_spn.value = getuserprop $PBRT_Helper "Render_Integ_IGI_Gather_spn" 
			
			PBRT_Render.Render_Integ_IrCach_MinWht_spn.value = getuserprop $PBRT_Helper "Render_Integ_IrCach_MinWht_spn" 
			PBRT_Render.Render_Integ_IrCach_MinSpac_spn.value = getuserprop $PBRT_Helper "Render_Integ_IrCach_MinSpac_spn" 
			PBRT_Render.Render_Integ_IrCach_MaxSpac_spn.value = getuserprop $PBRT_Helper "Render_Integ_IrCach_MaxSpac_spn" 
			PBRT_Render.Render_Integ_IrCach_MaxAngle_spn.value = getuserprop $PBRT_Helper "Render_Integ_IrCach_MaxAngle_spn" 
			PBRT_Render.Render_Integ_IrCach_SpecDepth_spn.value =getuserprop $PBRT_Helper "Render_Integ_IrCach_SpecDepth_spn"  
			PBRT_Render.Render_Integ_IrCach_indDepth_spn.value = getuserprop $PBRT_Helper "Render_Integ_IrCach_indDepth_spn" 
			PBRT_Render.Render_Integ_IrCach_Samples_spn.value = getuserprop $PBRT_Helper "Render_Integ_IrCach_Samples_spn" 
			
			PBRT_Render.Render_Integ_Path_Depth_spn.value = getuserprop $PBRT_Helper "Render_Integ_Path_Depth_spn" 
			
			PBRT_Render.Render_Integ_Photon_Caustic_spn.value = getuserprop $PBRT_Helper "Render_Integ_Photon_Caustic_spn" 
			PBRT_Render.Render_Integ_Photon_Indirect_spn.value = getuserprop $PBRT_Helper "Render_Integ_Photon_Indirect_spn" 
			PBRT_Render.Render_Integ_Photon_nUsed_spn.value = getuserprop $PBRT_Helper "Render_Integ_Photon_nUsed_spn" 
			PBRT_Render.Render_Integ_Photon_SpecDepth_spn.value = getuserprop $PBRT_Helper "Render_Integ_Photon_SpecDepth_spn" 
			PBRT_Render.Render_Integ_Photon_MaxPhoDepth_spn.value = getuserprop $PBRT_Helper "Render_Integ_Photon_MaxPhoDepth_spn" 
			PBRT_Render.Render_Integ_Photon_FinalGath_btn.state = getuserprop $PBRT_Helper "Render_Integ_Photon_FinalGath_btn" 
			PBRT_Render.Render_Integ_Photon_gatherSamp_spn.value = getuserprop $PBRT_Helper "Render_Integ_Photon_gatherSamp_spn" 
			PBRT_Render.Render_Integ_Photon_MaxDist_spn.value = getuserprop $PBRT_Helper "Render_Integ_Photon_MaxDist_spn" 
			PBRT_Render.Render_Integ_Photon_GathAng_spn.value = getuserprop $PBRT_Helper "Render_Integ_Photon_GathAng_spn" 
			
			PBRT_Render.Render_Integ_Single_Step_spn.value = getuserprop $PBRT_Helper "Render_Integ_Single_Step_spn" 
			
			PBRT_Render.Render_Integ_Whitted_MaDepth_spn.value = getuserprop $PBRT_Helper "Render_Integ_Whitted_MaDepth_spn" 
			
			PBRT_Render.Render_Vol_ddl.selection = getuserprop $PBRT_Helper "Render_Vol_ddl" 
			PBRT_Render.Render_Vol_Step_spn.value = getuserprop $PBRT_Helper "Render_Vol_Step_spn" 
			
			PBRT_Render.Render_Filter_ddl.selection = getuserprop $PBRT_Helper "Render_Filter_ddl" 
			
			PBRT_Render.Render_Filter_Box_X.value =getuserprop $PBRT_Helper "Render_Filter_Box_X"  
			PBRT_Render.Render_Filter_Box_Y.value = getuserprop $PBRT_Helper "Render_Filter_Box_Y" 
			
			PBRT_Render.Render_Filter_Gauss_X.value = getuserprop $PBRT_Helper "Render_Filter_Gauss_X" 
			PBRT_Render.Render_Filter_Gauss_Y.value = getuserprop $PBRT_Helper "Render_Filter_Gauss_Y" 
			PBRT_Render.Render_Filter_Gauss_Alpah.value = getuserprop $PBRT_Helper "Render_Filter_Gauss_Alpah" 
			
			PBRT_Render.Render_Filter_Mitch_X.value = getuserprop $PBRT_Helper "Render_Filter_Mitch_X" 
			PBRT_Render.Render_Filter_Mitch_Y.value = getuserprop $PBRT_Helper "Render_Filter_Mitch_Y" 
			PBRT_Render.Render_Filter_Mitch_Sharp.value = getuserprop $PBRT_Helper "Render_Filter_Mitch_Sharp" 
			
			PBRT_Render.Render_Filter_Sinc_X.value = getuserprop $PBRT_Helper "Render_Filter_Sinc_X" 
			PBRT_Render.Render_Filter_Sinc_Y.value = getuserprop $PBRT_Helper "Render_Filter_Sinc_Y" 
			PBRT_Render.Render_Filter_Sinc_Tau.value = getuserprop $PBRT_Helper "Render_Filter_Sinc_Tau" 
			
			PBRT_Render.Render_Filter_Tri_X.value = getuserprop $PBRT_Helper "Render_Filter_Tri_X" 
			PBRT_Render.Render_Filter_Tri_Y.value = getuserprop $PBRT_Helper "Render_Filter_Tri_Y" 

			PBRT_Render.exportObjs.state  = getuserprop $PBRT_Helper "exportObjs"  
			PBRT_Render.exportMaps.state  = getuserprop $PBRT_Helper "exportMaps" 
			PBRT_Render.envMult.value	 = getuserprop $PBRT_Helper "envMult" 
			PBRT_Render.envSamp.value = getuserprop $PBRT_Helper "envSamp" 
			
			PBRT_Render.useGlobalMat_chk.state = getuserprop $PBRT_Helper "useGlobalMat_chk" 	
		)
)

rollout PBRT_Render "Render" 
(
	dropdownList Render_Sampler_ddl "" 				pos:[27,23] 	width:127 	height:21 items:#("Low Discrepancy", "Random")
	spinner 		Render_Sampler_spn "Samples:" 	pos:[183,27] 	width:85 	height:16 range:[1,2048,4] type:#integer
	GroupBox 		Render_Sampler_Grp "Sampler" 	pos:[5,4] 		width:294 	height:43

	GroupBox 			Render_Rnd_grp 				"Render" 					pos:[5,51] 		width:294 	height:91
	dropdownList 	Render_Rnd_ddl 				"" 								pos:[15,70] 	width:98 	height:21 items:#("Aggregate Test", "Create Probes", "Metropolis", "Standard","Surface Points") selection:3

	spinner 		Render_Rnd_Aggre_Niters_spn "Niters:" pos:[103,101] width:85 height:16 range:[0,1e+006,100000] visible:false type:#integer
	
	spinner 		Render_Rnd_Probe_Lmax_spn 		"Lmax:" 						pos:[32,95] width:50 height:16 range:[0,100,4] visible:false type:#integer	
	spinner 		Render_Rnd_Probe_MaxDepth_spn 	"Max Depth:" 				pos:[51,120] width:50 height:16 range:[0,1000,5] visible:false type:#integer	
	spinner 		Render_Rnd_Probe_IndSamp_spn 	"Ind. Samples:" 			pos:[135,95] width:70 height:16 range:[0,100000,512] visible:false type:#integer
	spinner 		Render_Rnd_Probe_SampSpace_spn "Sam. Spac.:" 			pos:[155,120] width:60 height:16 range:[0,100000,1] visible:false type:#float	
	checkbutton 	Render_Rnd_Probe_Direct_btn 		"Direct Lighting:" 			pos:[208,63] width:85 height:17 visible:false toolTip:"Only Compute Indirect Illumination"	
	checkbutton 	Render_Rnd_Probe_Indirect_btn 		"Ind. Lighting:" 			pos:[208,84] width:85 height:17 visible:false toolTip:"Only Compute Indirect Illumination"		
	spinner 		Render_Rnd_Probe_time_spn 			"Time:" 						pos:[242,120] width:50 height:16 range:[0,100000,0] visible:false type:#float			
	
	spinner		 	Render_Rnd_Metro_Probab_spn 		"Prob:" 		pos:[30,97] width:45 height:16 range:[0,100,.25] visible:true type:#float
	spinner		 	Render_Rnd_Metro_Nsamples_spn 	"nSamp:" 		pos:[37,117] width:45 height:16 range:[0,100,0] visible:true type:#integer
	spinner		 	Render_Rnd_Metro_Samples_spn 	"Sam. Per Pixel:" 		pos:[126,97] width:85 height:16 range:[0,10000,100] visible:true type:#integer	
	spinner		 	Render_Rnd_Metro_BootStr_spn 		"Boot Str:" 		pos:[135,60] width:85 height:16 range:[0,100000,1000] visible:true type:#integer
	spinner 		Render_Rnd_Metro_DirectSamp 		"Direct Samp:" 		pos:[235,120] width:60 height:16 range:[1,1024,6] visible:true type:#integer toolTip:"Samples Per Pixel for separate direct lighting"			
	checkbutton 	Render_Rnd_Metro_DoDirect 			"Do Direct Sep:" pos:[222,85] width:72 height:17 checked:true visible:true toolTip:"Do Direct Lighitng Seperately"		
	spinner		 	Render_Rnd_Metro_MaxReje_spn 	"Max Reje:" 		pos:[155,80] width:65 height:16 range:[0,100000,512] visible:true type:#integer
	spinner 		Render_Rnd_Metro_MaxDepth 		"Max Depth:" 				pos:[130,117] width:60 height:16 range:[1,32,7] visible:true type:#integer	
	checkbutton 	Render_Rnd_Metro_Indirect_btn 		"Ind. Only:" 			pos:[225,60] width:70 height:17 visible:true toolTip:"Only Compute Indirect Illumination"
	
	spinner 		Render_Rnd_Surface_MinSam_spn 	"Min Sample Distance:" 		pos:[108,102] width:70 height:16 range:[0,100000,0.25] visible:false type:#float
		
	GroupBox 			Render_Inte_grp 									"Integrator" 				pos:[5,145] 	width:294 	height:143
	dropdownList 	Render_Integ_ddl 								"" 								pos:[71,163] 	width:159 	height:21 items:#("Ambient Occlusion", "Diffuse Seprt", \
							"Dipole Sub Surface", "Direct Lighting", "Emission", "Glossy Prt", "Instant Global Illumination", "Irradiance Cache", "Path Tracing", "Photon Map", \
						"Single", "Use Probes", "Whitted") selection:9
	
	spinner 			Render_Integ_AO_Samples_spn 				"Samples:" 					pos:[127,191] width:80 	height:16 range:[1,2048,32] visible:false type:#integer
	spinner 			Render_Integ_AO_MaxDist_spn 				"Max Distance:" 			pos:[92,212] 	width:115 	height:16 range:[1,1e+006,4] visible:false type:#float
	
	spinner 			Render_Integ_Diffuse_Max_spn 				"LMax:" 						pos:[105,191] width:91 	height:16 range:[1,2048,4] visible:false type:#integer
	spinner 			Render_Integ_Diffuse_Ns_spn 				"Ns:" 							pos:[117,212] width:80 	height:16 range:[1,4096,4] visible:false type:#integer
	
	spinner 			Render_Integ_Dipole_MaxDepth_spn 		"Max Depth:" 				pos:[104,192] width:91 	height:16 range:[1,2048,5] visible:false type:#integer
	spinner 			Render_Integ_Dipole_MaxError_spn 			"Max Error:" 				pos:[111,212] width:84 	height:16 range:[0.01,2048,0.05] visible:false type:#float
	spinner 			Render_Integ_Dipole_MinDist_spn 			"Min Distance:" 			pos:[94,231] 	width:101 	height:16 range:[0.01,2048,0.25] visible:false type:#float
	
	spinner 			Render_Integ_Direct_Depth_spn 			"Depth:" 					pos:[104,193] width:91 	height:16 range:[0,2048,8] visible:false	
	
	spinner 			Render_Integ_Emission_Step_spn 			"Step:" 						pos:[121,194] width:66 	height:16 range:[0.01,2048,1] visible:false type:#float
	
	spinner 			Render_Integ_Glossy_LMax_spn 				"LMax:" 						pos:[52,193] width:72 height:16 range:[1,2048,4] visible:false type:#integer
	spinner 			Render_Integ_Glossy_Ns_spn 				"Ns:" 							pos:[63,213] width:61 height:16 range:[1,2048,5] visible:false type:#integer
	checkbutton 		Render_Integ_Glossy_DotR_btn 				"Dot Transfer:" 			pos:[50,236] width:74 height:17  visible:false
	spinner 			Render_Integ_Glossy_Kd_spn 				"Kd:" 							pos:[191,195] width:61 height:16 range:[0.1,2048,0.5] visible:false type:#float
	spinner 			Render_Integ_Glossy_Ks_spn 					"Ks:" 							pos:[191,215] width:61 height:16 range:[0.1,2048,0.25] visible:false type:#float
	spinner 			Render_Integ_Glossy_Rough_spn 			"Roughness:" 				pos:[162,236] width:90 height:16 range:[0.1,2048,10] visible:false type:#float
	
	spinner 			Render_Integ_IGI_LightPath_spn 			"Light Paths:" 				pos:[50,189] width:90 height:16 range:[1,2048,64] visible:false type:#integer
	spinner 			Render_Integ_IGI_LightSet_spn 				"Light Sets:" 				pos:[50,209] width:90 height:16 range:[1,2048,4] visible:false type:#integer
	spinner 			Render_Integ_IGI_MinDist_spn 				"Min Distance:" 			pos:[50,229] width:90 height:16 range:[1,2048,0.1] visible:false type:#float
	spinner 			Render_Integ_IGI_RRThrs_spn 				"RR Threshold:" 			pos:[50,248] width:90 height:16 range:[1,2048,0.0001] visible:false type:#float
	spinner 			Render_Integ_IGI_MaxDepth_spn 			"Max Depth:" 				pos:[191,190] width:85 height:16 range:[1,2048,5] visible:false type:#integer
	spinner 			Render_Integ_IGI_GLimit_spn 				"G Limit:" 					pos:[201,210] width:76 height:16 range:[0.1,2048,10] visible:false type:#float
	spinner 			Render_Integ_IGI_Gather_spn 				"Gather Samples:" 		pos:[192,229] width:85 height:16 range:[1,2048,16] visible:false type:#integer
	
	spinner 			Render_Integ_IrCach_MinWht_spn 			"Min Weight:" 				pos:[45,192] width:90 height:16 range:[0.1,2048,0.5] visible:false type:#float
	spinner 			Render_Integ_IrCach_MinSpac_spn 		"Min Spacing:" 			pos:[45,212] width:90 height:16 range:[0.1,2048,2.5] visible:false type:#float
	spinner 			Render_Integ_IrCach_MaxSpac_spn 		"Max Spacing:" 			pos:[45,232] width:90 height:16 range:[0.1,2048,15] visible:false type:#float
	spinner 			Render_Integ_IrCach_MaxAngle_spn 		"Max Angle:" 				pos:[215,190] width:75 height:16 range:[0.1,2048,10] visible:false type:#float
	spinner 			Render_Integ_IrCach_SpecDepth_spn 		"Max Specular Depth:" 	pos:[206,210] width:85 height:16 range:[0.1,2048,5] visible:false type:#integer
	spinner 			Render_Integ_IrCach_indDepth_spn 		"Max Indirect Depth:" 	pos:[206,230] width:85 height:16 range:[0.1,2048,3] visible:false type:#integer
	spinner 			Render_Integ_IrCach_Samples_spn 			"Samples:" 					pos:[2120,250] width:70 height:16 range:[0.1,2048,64] visible:false type:#integer
	
	spinner 			Render_Integ_Path_Depth_spn 				"Depth:" 					pos:[110,193] width:85 height:16 range:[1,2048,5] visible:true type:#integer
	
	spinner 			Render_Integ_Photon_Caustic_spn 			"Caustic:" 					pos:[50,190] width:99 height:16 range:[1,5e+006,20000] visible:false type:#integer
	spinner 			Render_Integ_Photon_Indirect_spn 			"Indirect:" 					pos:[50,209] width:99 height:16 range:[1,5e+006,100000] visible:false type:#integer
	spinner 			Render_Integ_Photon_nUsed_spn 			"Nused:" 					pos:[70,227] width:80 height:16 range:[0,100000,50] visible:false type:#integer	
	spinner 			Render_Integ_Photon_SpecDepth_spn 	"Max Specular Depth:" 	pos:[80,246] width:70 height:16 range:[1,1024,5] visible:false type:#integer
	spinner 			Render_Integ_Photon_MaxPhoDepth_spn "Max Photon Depth:" 		pos:[81,265] width:68 height:16 range:[1,1024,10] visible:false type:#integer	
	checkbutton		Render_Integ_Photon_FinalGath_btn 		"Final Gather:" 			pos:[215,189] width:70 height:17  visible:false
	spinner 			Render_Integ_Photon_gatherSamp_spn 	"Gather Samples:" 		pos:[205,210] width:90 height:16 range:[1,1024,32] visible:false type:#integer
	spinner 			Render_Integ_Photon_MaxDist_spn 			"Max Distance:" 			pos:[205,231] width:90 height:16 range:[0.01,1024,0.1] visible:false type:#float
	spinner 			Render_Integ_Photon_GathAng_spn 		"Gather Angle:" 			pos:[206,251] width:90 height:16 range:[0.01,1024,10] visible:false type:#float
	
	spinner 			Render_Integ_Single_Step_spn 				"Step Size:" 				pos:[108,195] width:86 height:16 range:[0.01,1024,0.1] visible:false type:#float
	
	spinner 			Render_Integ_Whitted_MaDepth_spn 		"Max Depth:" 				pos:[110,193] width:86 height:16 range:[1,1024,5] visible:false type:#integer
	
	GroupBox 			Render_Vol_grp 			"Volume Int" 	pos:[5,292] width:294 height:49
	dropdownList 	Render_Vol_ddl 			"" 					pos:[23,311] width:127 height:21 items:#("Single", "Emission")
	spinner 			Render_Vol_Step_spn 	"Step Size:" 	pos:[176,314] width:103 height:16 range:[0,100,1] type:#float
	
	GroupBox 			Render_Filter_grp 			"Filter" 			pos:[5,343] width:294 height:95
	dropdownList 	Render_Filter_ddl 			"" 					pos:[82,358] width:127 height:21 items:#("box", "gaussian", "mitchell", "sinc", "triangle") selection:3
	
	spinner 			Render_Filter_Box_X 		"Width X:" 		pos:[46,390] width:90 height:16 range:[0,10,0.5] type:#float visible:false
	spinner 			Render_Filter_Box_Y 		"Width Y:" 		pos:[163,391] width:90 height:16 range:[0,10,0.5] type:#float visible:false
	
	spinner 			Render_Filter_Gauss_X 		"Width X:" 		pos:[46,391] width:90 height:16 range:[0,10,0.5] type:#float visible:false
	spinner 			Render_Filter_Gauss_Y 		"Width Y:" 		pos:[162,392] width:90 height:16 range:[0,10,0.5] type:#float visible:false
	spinner 			Render_Filter_Gauss_Alpah "Alpah:" 		pos:[65,413] width:71 height:16 range:[0,10,2] type:#float visible:false
	
	spinner 			Render_Filter_Mitch_X 		"Width X:" 		pos:[46,391] width:90 height:16 range:[0,10,0.5] type:#float
	spinner 			Render_Filter_Mitch_Y 		"Width Y:" 		pos:[162,391] width:90 height:16 range:[0,10,0.5] type:#float
	spinner 			Render_Filter_Mitch_Sharp "Sharpness:" 	pos:[45,415] width:90 height:16 range:[0,1,0.33] type:#float
	
	spinner 			Render_Filter_Sinc_X 		"Width X:" 		pos:[46,391] width:90 height:16 range:[0,10,0.5] type:#float visible:false
	spinner 			Render_Filter_Sinc_Y 		"Width Y:" 		pos:[162,391] width:90 height:16 range:[0,10,0.5] type:#float visible:false
	spinner 			Render_Filter_Sinc_Tau 	"Tau:" 			pos:[77,416] width:57 height:16 range:[0,1,3] type:#float visible:false
	
	spinner 			Render_Filter_Tri_X 			"Width X:" 		pos:[46,391] width:90 height:16 range:[0,10,0.5] type:#float visible:false
	spinner 			Render_Filter_Tri_Y 			"Width Y:" 		pos:[162,391] width:90 height:16 range:[0,10,0.5] type:#float visible:false
	
	
	button RenderIT	"RENDER" pos:[5, 440] width:140 height:17 enabled:false
	button ExportIT	"EXPORT" pos:[152, 440] width:140 height:17
	checkbox exportObjs "Export Geometry" pos:[25, 460] width:100 height:17 checked:true  
	checkbox exportMaps "Export Materials" pos:[150, 460] width:100 height:17 checked:true  
	spinner envMult "Env.Multiplier" pos:[45, 480] range:[1.0, 100.0, 1.0] width:100 type:#float	
	spinner envSamp "Env. Samples" pos:[200, 480] range:[1, 200, 32] width:90 type:#integer
	
	checkbox 				useGlobalMat_chk "Override Global Mat." pos:[5, 500] width:120 height:17 checked:false
	materialbutton 		GloablMaterial_btn	"Pick Material" pos:[130, 500] width:170 height:17 tooltip:"Select Override Map from the Material Editor"
	
	progressBar export_prog "" pos:[5,520] width:290 height:7 color:(color 243 103 0) 
	
	on PBRT_Render open do PBRT_GetProperties()
	on PBRT_Render close do PBRT_SaveProperties()
		
		on exportMaps changed state do (
	--		doneMats = #() -- cler the array so it can be recreated
		)
	
	on GloablMaterial_btn picked val do (
		GlobalMaterial = val
		GloablMaterial_btn.text = GlobalMaterial .name
	)
	
	on RenderIT pressed do (
				
--		if exportMaps.checked then doneMats = #() -- cler the array so it can be recreated
		doneMats = #() -- array to hold the materials that have been already processed
		TextArray = #()	-- array to hold the proccesed textures
		
	--	Local dirPath = PBRT_System.PBRT_exe_Path_edt.text 
		Local theFileOnly = filterstring PBRT_System.PBRT_SavePath_edt.text "\\"		
		
		-- we need to check if the user provided file extension so we dont duplicat it
		Local theName = theFileOnly[theFileOnly.count] 	-- name pat of the complete string with or without extension, previously we already spli it by \
		Local remExt =  filterstring theName "."	-- split it again now by .

		-- if it has cont more then 1 it means it has an extension
		if remExt.count > 1 then theFileOnly = remExt[1] -- assign the element 1, removing the extension
			else theFileOnly = remExt[remExt.count] 		-- the user did not provide the extension
	
		Local dirPath = PBRT_System.PBRT_SavePath_edt.text 
		Only = filterstring dirPath "\\"
		newTxt = substring dirPath 1 (dirPath.count - only[only.count].count - 1)
		dirPath = newTxt + "\\"
		
		setWaitCursor()
		Pbrt_Export()

		thefile = PBRT_System.PBRT_SavePath_edt.text + ".pbrt"
--		theCommand = "\"\"" + dirPath + "pbrt.exe\" \"" + thefile + "\"\"" --"\"\"" + dirPath + "pbrt.exe\" \"" + (theFileOnly[theFileOnly.count]  + ".pbrt\"\"")
--		theCommand = "\"\"" + dirPath + "pbrt.exe\" \"" + (theFileOnly[theFileOnly.count]  + ".pbrt\"\"")
		
		theCommand = "\"\"" + dirPath + "pbrt.exe\" \"" + (theFileOnly + ".pbrt\"\"")

		HiddenDOSCommand theCommand startpath:dirPath --donotwait:true

		setArrowCursor()	

		
		Local out =  Viewer + "exrdisplay.exe" +  "\"" 
print ("\"" +out  + (" \"" + dirpath + theFileOnly  + ".exr\"") + "\"")
		DOSCommand  ("\"" +out  + (" \"" + dirpath + theFileOnly   + ".exr\"") + "\"")				
	)

	on ExportIT pressed do (

--		if exportMaps.checked then doneMats = #() -- cler the array so it can be recreated
		doneMats = #() -- array to hold the materials that have been already processed
		TextArray = #()	-- array to hold the proccesed textures
		
		setWaitCursor()
		Pbrt_Export()
		setArrowCursor()	
		
	)
	
	on Render_Rnd_ddl  selected arg do (
		case Render_Rnd_ddl.selection of (
			1: -- Aggregate Test
			(
			Render_Rnd_Aggre_Niters_spn.visible = true
			Render_Rnd_Probe_Direct_btn.visible = false
			Render_Rnd_Probe_Indirect_btn.visible = false	
			Render_Rnd_Probe_Lmax_spn.visible = false
			Render_Rnd_Probe_IndSamp_spn.visible = false
			Render_Rnd_Probe_MaxDepth_spn.visible = false
			Render_Rnd_Probe_SampSpace_spn.visible = false
			Render_Rnd_Probe_time_spn.visible = false					
			Render_Rnd_Metro_Probab_spn.visible = false
			Render_Rnd_Metro_Nsamples_spn.visible = false
			Render_Rnd_Metro_Samples_spn.visible = false
			Render_Rnd_Metro_BootStr_spn.visible = false
			Render_Rnd_Metro_DirectSamp.visible = false		
			Render_Rnd_Metro_DoDirect.visible = false	
			Render_Rnd_Metro_MaxReje_spn.visible = false
			Render_Rnd_Metro_MaxDepth.visible = false	
			Render_Rnd_Metro_Indirect_btn.visible = false				
			Render_Rnd_Surface_MinSam_spn.visible = false
			)
			2: -- Create Probes
			(
			Render_Rnd_Aggre_Niters_spn.visible = false
			Render_Rnd_Probe_Direct_btn.visible = true
			Render_Rnd_Probe_Indirect_btn.visible = true	
			Render_Rnd_Probe_Lmax_spn.visible = true
			Render_Rnd_Probe_IndSamp_spn.visible = true
			Render_Rnd_Probe_MaxDepth_spn.visible = true
			Render_Rnd_Probe_SampSpace_spn.visible = true
			Render_Rnd_Probe_time_spn.visible = true				
			Render_Rnd_Metro_Probab_spn.visible = false
			Render_Rnd_Metro_Nsamples_spn.visible = false
			Render_Rnd_Metro_Samples_spn.visible = false
			Render_Rnd_Metro_BootStr_spn.visible = false
			Render_Rnd_Metro_DirectSamp.visible = false		
			Render_Rnd_Metro_DoDirect.visible = false	
			Render_Rnd_Metro_MaxReje_spn.visible = false
			Render_Rnd_Metro_MaxDepth.visible = false	
			Render_Rnd_Metro_Indirect_btn.visible = false
			Render_Rnd_Surface_MinSam_spn.visible = false				
			)
			3: -- Metropolis
			(
			Render_Rnd_Aggre_Niters_spn.visible = false
			Render_Rnd_Probe_Direct_btn.visible = false
			Render_Rnd_Probe_Indirect_btn.visible = false	
			Render_Rnd_Probe_Lmax_spn.visible = false
			Render_Rnd_Probe_IndSamp_spn.visible = false
			Render_Rnd_Probe_MaxDepth_spn.visible = false
			Render_Rnd_Probe_SampSpace_spn.visible = false
			Render_Rnd_Probe_time_spn.visible = false				
			Render_Rnd_Metro_Probab_spn.visible = true
			Render_Rnd_Metro_Nsamples_spn.visible = true
			Render_Rnd_Metro_Samples_spn.visible = true
			Render_Rnd_Metro_BootStr_spn.visible = true
			Render_Rnd_Metro_DirectSamp.visible = true		
			Render_Rnd_Metro_DoDirect.visible = true	
			Render_Rnd_Metro_MaxReje_spn.visible = true
			Render_Rnd_Metro_MaxDepth.visible = true
			Render_Rnd_Metro_Indirect_btn.visible = true
			Render_Rnd_Surface_MinSam_spn.visible = false					
			)
			4: -- Standard
			(
			Render_Rnd_Aggre_Niters_spn.visible = false
			Render_Rnd_Probe_Direct_btn.visible = false
			Render_Rnd_Probe_Indirect_btn.visible = false	
			Render_Rnd_Probe_Lmax_spn.visible = false
			Render_Rnd_Probe_IndSamp_spn.visible = false
			Render_Rnd_Probe_MaxDepth_spn.visible = false
			Render_Rnd_Probe_SampSpace_spn.visible = false
			Render_Rnd_Probe_time_spn.visible = false			
			Render_Rnd_Metro_Probab_spn.visible = false
			Render_Rnd_Metro_Nsamples_spn.visible = false
			Render_Rnd_Metro_Samples_spn.visible = false
			Render_Rnd_Metro_BootStr_spn.visible = false
			Render_Rnd_Metro_DirectSamp.visible = false		
			Render_Rnd_Metro_DoDirect.visible = false	
			Render_Rnd_Metro_MaxReje_spn.visible = false
			Render_Rnd_Metro_MaxDepth.visible = false	
			Render_Rnd_Metro_Indirect_btn.visible = false
			Render_Rnd_Surface_MinSam_spn.visible = false					
			)
			5: -- Surface Points
			(
			Render_Rnd_Aggre_Niters_spn.visible = false
			Render_Rnd_Probe_Direct_btn.visible = false
			Render_Rnd_Probe_Indirect_btn.visible = false	
			Render_Rnd_Probe_Lmax_spn.visible = false
			Render_Rnd_Probe_IndSamp_spn.visible = false
			Render_Rnd_Probe_MaxDepth_spn.visible = false
			Render_Rnd_Probe_SampSpace_spn.visible = false
			Render_Rnd_Probe_time_spn.visible = false				
			Render_Rnd_Metro_Probab_spn.visible = false
			Render_Rnd_Metro_Nsamples_spn.visible = false
			Render_Rnd_Metro_Samples_spn.visible = false
			Render_Rnd_Metro_BootStr_spn.visible = false
			Render_Rnd_Metro_DirectSamp.visible = false		
			Render_Rnd_Metro_DoDirect.visible = false	
			Render_Rnd_Metro_MaxReje_spn.visible = false
			Render_Rnd_Metro_MaxDepth.visible = false	
			Render_Rnd_Metro_Indirect_btn.visible = false
			Render_Rnd_Surface_MinSam_spn.visible = true					
			)
		)
	)
	
	on Render_Integ_ddl selected arg do (
		case Render_Integ_ddl.selection of (
			1:	-- Ambient Occlusion
			(
			Render_Integ_AO_Samples_spn.visible = true
			Render_Integ_AO_MaxDist_spn.visible = true				
			Render_Integ_Diffuse_Max_spn.visible = false
			Render_Integ_Diffuse_Ns_spn.visible = false				
			Render_Integ_Dipole_MaxDepth_spn.visible = false
			Render_Integ_Dipole_MaxError_spn.visible = false
			Render_Integ_Dipole_MinDist_spn.visible = false
			Render_Integ_Direct_Depth_spn.visible = false					
			Render_Integ_Emission_Step_spn.visible = false				
			Render_Integ_Glossy_LMax_spn.visible = false
			Render_Integ_Glossy_Ns_spn.visible = false
			Render_Integ_Glossy_DotR_btn.visible = false
			Render_Integ_Glossy_Kd_spn.visible = false
			Render_Integ_Glossy_Ks_spn.visible = false
			Render_Integ_Glossy_Rough_spn.visible = false				
			Render_Integ_IGI_LightPath_spn.visible = false
			Render_Integ_IGI_LightSet_spn.visible = false
			Render_Integ_IGI_MinDist_spn.visible = false
			Render_Integ_IGI_RRThrs_spn.visible = false
			Render_Integ_IGI_MaxDepth_spn.visible = false
			Render_Integ_IGI_GLimit_spn.visible = false
			Render_Integ_IGI_Gather_spn.visible = false			
			Render_Integ_IrCach_MinWht_spn.visible = false
			Render_Integ_IrCach_MinSpac_spn.visible = false
			Render_Integ_IrCach_MaxSpac_spn.visible = false
			Render_Integ_IrCach_MaxAngle_spn.visible = false
			Render_Integ_IrCach_SpecDepth_spn.visible = false
			Render_Integ_IrCach_indDepth_spn.visible = false
			Render_Integ_IrCach_Samples_spn.visible = false				
			Render_Integ_Path_Depth_spn.visible = false				
			Render_Integ_Photon_Caustic_spn.visible = false
			Render_Integ_Photon_Indirect_spn.visible = false
			Render_Integ_Photon_nUsed_spn.visible = false
			Render_Integ_Photon_SpecDepth_spn.visible = false
			Render_Integ_Photon_MaxPhoDepth_spn.visible = false
			Render_Integ_Photon_FinalGath_btn.visible = false
			Render_Integ_Photon_gatherSamp_spn.visible = false
			Render_Integ_Photon_MaxDist_spn.visible = false
			Render_Integ_Photon_GathAng_spn.visible = false
			Render_Integ_Single_Step_spn.visible = false				
			Render_Integ_Whitted_MaDepth_spn.visible = false			
			)
			2: -- Diffuse Seprt
			(
			Render_Integ_AO_Samples_spn.visible = false
			Render_Integ_AO_MaxDist_spn.visible = false				
			Render_Integ_Diffuse_Max_spn.visible = true
			Render_Integ_Diffuse_Ns_spn.visible = true				
			Render_Integ_Dipole_MaxDepth_spn.visible = false
			Render_Integ_Dipole_MaxError_spn.visible = false
			Render_Integ_Dipole_MinDist_spn.visible = false
			Render_Integ_Direct_Depth_spn.visible = false					
			Render_Integ_Emission_Step_spn.visible = false				
			Render_Integ_Glossy_LMax_spn.visible = false
			Render_Integ_Glossy_Ns_spn.visible = false
			Render_Integ_Glossy_DotR_btn.visible = false
			Render_Integ_Glossy_Kd_spn.visible = false
			Render_Integ_Glossy_Ks_spn.visible = false
			Render_Integ_Glossy_Rough_spn.visible = false				
			Render_Integ_IGI_LightPath_spn.visible = false
			Render_Integ_IGI_LightSet_spn.visible = false
			Render_Integ_IGI_MinDist_spn.visible = false
			Render_Integ_IGI_RRThrs_spn.visible = false
			Render_Integ_IGI_MaxDepth_spn.visible = false
			Render_Integ_IGI_GLimit_spn.visible = false
			Render_Integ_IGI_Gather_spn.visible = false			
			Render_Integ_IrCach_MinWht_spn.visible = false
			Render_Integ_IrCach_MinSpac_spn.visible = false
			Render_Integ_IrCach_MaxSpac_spn.visible = false
			Render_Integ_IrCach_MaxAngle_spn.visible = false
			Render_Integ_IrCach_SpecDepth_spn.visible = false
			Render_Integ_IrCach_indDepth_spn.visible = false
			Render_Integ_IrCach_Samples_spn.visible = false				
			Render_Integ_Path_Depth_spn.visible = false				
			Render_Integ_Photon_Caustic_spn.visible = false
			Render_Integ_Photon_Indirect_spn.visible = false
			Render_Integ_Photon_nUsed_spn.visible = false
			Render_Integ_Photon_SpecDepth_spn.visible = false
			Render_Integ_Photon_MaxPhoDepth_spn.visible = false
			Render_Integ_Photon_FinalGath_btn.visible = false
			Render_Integ_Photon_gatherSamp_spn.visible = false
			Render_Integ_Photon_MaxDist_spn.visible = false
			Render_Integ_Photon_GathAng_spn.visible = false
			Render_Integ_Single_Step_spn.visible = false				
			Render_Integ_Whitted_MaDepth_spn.visible = false					
			)
			3: -- Dipole Sub Surface
			(
			Render_Integ_AO_Samples_spn.visible = false
			Render_Integ_AO_MaxDist_spn.visible = false				
			Render_Integ_Diffuse_Max_spn.visible = false
			Render_Integ_Diffuse_Ns_spn.visible = false				
			Render_Integ_Dipole_MaxDepth_spn.visible = true
			Render_Integ_Dipole_MaxError_spn.visible = true
			Render_Integ_Dipole_MinDist_spn.visible = true	
			Render_Integ_Direct_Depth_spn.visible = false					
			Render_Integ_Emission_Step_spn.visible = false				
			Render_Integ_Glossy_LMax_spn.visible = false
			Render_Integ_Glossy_Ns_spn.visible = false
			Render_Integ_Glossy_DotR_btn.visible = false
			Render_Integ_Glossy_Kd_spn.visible = false
			Render_Integ_Glossy_Ks_spn.visible = false
			Render_Integ_Glossy_Rough_spn.visible = false				
			Render_Integ_IGI_LightPath_spn.visible = false
			Render_Integ_IGI_LightSet_spn.visible = false
			Render_Integ_IGI_MinDist_spn.visible = false
			Render_Integ_IGI_RRThrs_spn.visible = false
			Render_Integ_IGI_MaxDepth_spn.visible = false
			Render_Integ_IGI_GLimit_spn.visible = false
			Render_Integ_IGI_Gather_spn.visible = false			
			Render_Integ_IrCach_MinWht_spn.visible = false
			Render_Integ_IrCach_MinSpac_spn.visible = false
			Render_Integ_IrCach_MaxSpac_spn.visible = false
			Render_Integ_IrCach_MaxAngle_spn.visible = false
			Render_Integ_IrCach_SpecDepth_spn.visible = false
			Render_Integ_IrCach_indDepth_spn.visible = false
			Render_Integ_IrCach_Samples_spn.visible = false				
			Render_Integ_Path_Depth_spn.visible = false				
			Render_Integ_Photon_Caustic_spn.visible = false
			Render_Integ_Photon_Indirect_spn.visible = false
			Render_Integ_Photon_nUsed_spn.visible = false
			Render_Integ_Photon_SpecDepth_spn.visible = false
			Render_Integ_Photon_MaxPhoDepth_spn.visible = false
			Render_Integ_Photon_FinalGath_btn.visible = false
			Render_Integ_Photon_gatherSamp_spn.visible = false
			Render_Integ_Photon_MaxDist_spn.visible = false
			Render_Integ_Photon_GathAng_spn.visible = false
			Render_Integ_Single_Step_spn.visible = false				
			Render_Integ_Whitted_MaDepth_spn.visible = false				
			)
			4: -- Direct Lighting
			(
			Render_Integ_AO_Samples_spn.visible = false
			Render_Integ_AO_MaxDist_spn.visible = false				
			Render_Integ_Diffuse_Max_spn.visible = false
			Render_Integ_Diffuse_Ns_spn.visible = false				
			Render_Integ_Dipole_MaxDepth_spn.visible = false	
			Render_Integ_Dipole_MaxError_spn.visible = false	
			Render_Integ_Dipole_MinDist_spn.visible = false	
			Render_Integ_Direct_Depth_spn.visible = true					
			Render_Integ_Emission_Step_spn.visible = false				
			Render_Integ_Glossy_LMax_spn.visible = false
			Render_Integ_Glossy_Ns_spn.visible = false
			Render_Integ_Glossy_DotR_btn.visible = false
			Render_Integ_Glossy_Kd_spn.visible = false
			Render_Integ_Glossy_Ks_spn.visible = false
			Render_Integ_Glossy_Rough_spn.visible = false				
			Render_Integ_IGI_LightPath_spn.visible = false
			Render_Integ_IGI_LightSet_spn.visible = false
			Render_Integ_IGI_MinDist_spn.visible = false
			Render_Integ_IGI_RRThrs_spn.visible = false
			Render_Integ_IGI_MaxDepth_spn.visible = false
			Render_Integ_IGI_GLimit_spn.visible = false
			Render_Integ_IGI_Gather_spn.visible = false			
			Render_Integ_IrCach_MinWht_spn.visible = false
			Render_Integ_IrCach_MinSpac_spn.visible = false
			Render_Integ_IrCach_MaxSpac_spn.visible = false
			Render_Integ_IrCach_MaxAngle_spn.visible = false
			Render_Integ_IrCach_SpecDepth_spn.visible = false
			Render_Integ_IrCach_indDepth_spn.visible = false
			Render_Integ_IrCach_Samples_spn.visible = false				
			Render_Integ_Path_Depth_spn.visible = false				
			Render_Integ_Photon_Caustic_spn.visible = false
			Render_Integ_Photon_Indirect_spn.visible = false
			Render_Integ_Photon_nUsed_spn.visible = false
			Render_Integ_Photon_SpecDepth_spn.visible = false
			Render_Integ_Photon_MaxPhoDepth_spn.visible = false
			Render_Integ_Photon_FinalGath_btn.visible = false
			Render_Integ_Photon_gatherSamp_spn.visible = false
			Render_Integ_Photon_MaxDist_spn.visible = false
			Render_Integ_Photon_GathAng_spn.visible = false
			Render_Integ_Single_Step_spn.visible = false				
			Render_Integ_Whitted_MaDepth_spn.visible = false				
			)
			5: -- Emission
			(
			Render_Integ_AO_Samples_spn.visible = false
			Render_Integ_AO_MaxDist_spn.visible = false				
			Render_Integ_Diffuse_Max_spn.visible = false
			Render_Integ_Diffuse_Ns_spn.visible = false				
			Render_Integ_Dipole_MaxDepth_spn.visible = false	
			Render_Integ_Dipole_MaxError_spn.visible = false	
			Render_Integ_Dipole_MinDist_spn.visible = false	
			Render_Integ_Direct_Depth_spn.visible = false					
			Render_Integ_Emission_Step_spn.visible = true				
			Render_Integ_Glossy_LMax_spn.visible = false
			Render_Integ_Glossy_Ns_spn.visible = false
			Render_Integ_Glossy_DotR_btn.visible = false
			Render_Integ_Glossy_Kd_spn.visible = false
			Render_Integ_Glossy_Ks_spn.visible = false
			Render_Integ_Glossy_Rough_spn.visible = false				
			Render_Integ_IGI_LightPath_spn.visible = false
			Render_Integ_IGI_LightSet_spn.visible = false
			Render_Integ_IGI_MinDist_spn.visible = false
			Render_Integ_IGI_RRThrs_spn.visible = false
			Render_Integ_IGI_MaxDepth_spn.visible = false
			Render_Integ_IGI_GLimit_spn.visible = false
			Render_Integ_IGI_Gather_spn.visible = false			
			Render_Integ_IrCach_MinWht_spn.visible = false
			Render_Integ_IrCach_MinSpac_spn.visible = false
			Render_Integ_IrCach_MaxSpac_spn.visible = false
			Render_Integ_IrCach_MaxAngle_spn.visible = false
			Render_Integ_IrCach_SpecDepth_spn.visible = false
			Render_Integ_IrCach_indDepth_spn.visible = false
			Render_Integ_IrCach_Samples_spn.visible = false				
			Render_Integ_Path_Depth_spn.visible = false				
			Render_Integ_Photon_Caustic_spn.visible = false
			Render_Integ_Photon_Indirect_spn.visible = false
			Render_Integ_Photon_nUsed_spn.visible = false
			Render_Integ_Photon_SpecDepth_spn.visible = false
			Render_Integ_Photon_MaxPhoDepth_spn.visible = false
			Render_Integ_Photon_FinalGath_btn.visible = false
			Render_Integ_Photon_gatherSamp_spn.visible = false
			Render_Integ_Photon_MaxDist_spn.visible = false
			Render_Integ_Photon_GathAng_spn.visible = false
			Render_Integ_Single_Step_spn.visible = false				
			Render_Integ_Whitted_MaDepth_spn.visible = false					
			)
			6: -- Glossy Prt
			(
			Render_Integ_AO_Samples_spn.visible = false
			Render_Integ_AO_MaxDist_spn.visible = false				
			Render_Integ_Diffuse_Max_spn.visible = false
			Render_Integ_Diffuse_Ns_spn.visible = false				
			Render_Integ_Dipole_MaxDepth_spn.visible = false	
			Render_Integ_Dipole_MaxError_spn.visible = false	
			Render_Integ_Dipole_MinDist_spn.visible = false	
			Render_Integ_Direct_Depth_spn.visible = false					
			Render_Integ_Emission_Step_spn.visible = false				
			Render_Integ_Glossy_LMax_spn.visible = true
			Render_Integ_Glossy_Ns_spn.visible = true
			Render_Integ_Glossy_DotR_btn.visible = true
			Render_Integ_Glossy_Kd_spn.visible = true
			Render_Integ_Glossy_Ks_spn.visible = true
			Render_Integ_Glossy_Rough_spn.visible = true				
			Render_Integ_IGI_LightPath_spn.visible = false
			Render_Integ_IGI_LightSet_spn.visible = false
			Render_Integ_IGI_MinDist_spn.visible = false
			Render_Integ_IGI_RRThrs_spn.visible = false
			Render_Integ_IGI_MaxDepth_spn.visible = false
			Render_Integ_IGI_GLimit_spn.visible = false
			Render_Integ_IGI_Gather_spn.visible = false			
			Render_Integ_IrCach_MinWht_spn.visible = false
			Render_Integ_IrCach_MinSpac_spn.visible = false
			Render_Integ_IrCach_MaxSpac_spn.visible = false
			Render_Integ_IrCach_MaxAngle_spn.visible = false
			Render_Integ_IrCach_SpecDepth_spn.visible = false
			Render_Integ_IrCach_indDepth_spn.visible = false
			Render_Integ_IrCach_Samples_spn.visible = false				
			Render_Integ_Path_Depth_spn.visible = false				
			Render_Integ_Photon_Caustic_spn.visible = false
			Render_Integ_Photon_Indirect_spn.visible = false
			Render_Integ_Photon_nUsed_spn.visible = false
			Render_Integ_Photon_SpecDepth_spn.visible = false
			Render_Integ_Photon_MaxPhoDepth_spn.visible = false
			Render_Integ_Photon_FinalGath_btn.visible = false
			Render_Integ_Photon_gatherSamp_spn.visible = false
			Render_Integ_Photon_MaxDist_spn.visible = false
			Render_Integ_Photon_GathAng_spn.visible = false
			Render_Integ_Single_Step_spn.visible = false				
			Render_Integ_Whitted_MaDepth_spn.visible = false					
			)
			7: -- Instant Global Illumination
			(
			Render_Integ_AO_Samples_spn.visible = false
			Render_Integ_AO_MaxDist_spn.visible = false				
			Render_Integ_Diffuse_Max_spn.visible = false
			Render_Integ_Diffuse_Ns_spn.visible = false				
			Render_Integ_Dipole_MaxDepth_spn.visible = false	
			Render_Integ_Dipole_MaxError_spn.visible = false	
			Render_Integ_Dipole_MinDist_spn.visible = false	
			Render_Integ_Direct_Depth_spn.visible = false					
			Render_Integ_Emission_Step_spn.visible = false				
			Render_Integ_Glossy_LMax_spn.visible = false
			Render_Integ_Glossy_Ns_spn.visible = false
			Render_Integ_Glossy_DotR_btn.visible = false
			Render_Integ_Glossy_Kd_spn.visible = false
			Render_Integ_Glossy_Ks_spn.visible = false
			Render_Integ_Glossy_Rough_spn.visible = false				
			Render_Integ_IGI_LightPath_spn.visible = true
			Render_Integ_IGI_LightSet_spn.visible = true
			Render_Integ_IGI_MinDist_spn.visible = true
			Render_Integ_IGI_RRThrs_spn.visible = true
			Render_Integ_IGI_MaxDepth_spn.visible = true
			Render_Integ_IGI_GLimit_spn.visible = true
			Render_Integ_IGI_Gather_spn.visible = true		
			Render_Integ_IrCach_MinWht_spn.visible = false
			Render_Integ_IrCach_MinSpac_spn.visible = false
			Render_Integ_IrCach_MaxSpac_spn.visible = false
			Render_Integ_IrCach_MaxAngle_spn.visible = false
			Render_Integ_IrCach_SpecDepth_spn.visible = false
			Render_Integ_IrCach_indDepth_spn.visible = false
			Render_Integ_IrCach_Samples_spn.visible = false				
			Render_Integ_Path_Depth_spn.visible = false				
			Render_Integ_Photon_Caustic_spn.visible = false
			Render_Integ_Photon_Indirect_spn.visible = false
			Render_Integ_Photon_nUsed_spn.visible = false
			Render_Integ_Photon_SpecDepth_spn.visible = false
			Render_Integ_Photon_MaxPhoDepth_spn.visible = false
			Render_Integ_Photon_FinalGath_btn.visible = false
			Render_Integ_Photon_gatherSamp_spn.visible = false
			Render_Integ_Photon_MaxDist_spn.visible = false
			Render_Integ_Photon_GathAng_spn.visible = false
			Render_Integ_Single_Step_spn.visible = false				
			Render_Integ_Whitted_MaDepth_spn.visible = false					
			)
			8: -- Irradiance Cache
			(
			Render_Integ_AO_Samples_spn.visible = false
			Render_Integ_AO_MaxDist_spn.visible = false				
			Render_Integ_Diffuse_Max_spn.visible = false
			Render_Integ_Diffuse_Ns_spn.visible = false				
			Render_Integ_Dipole_MaxDepth_spn.visible = false	
			Render_Integ_Dipole_MaxError_spn.visible = false	
			Render_Integ_Dipole_MinDist_spn.visible = false	
			Render_Integ_Direct_Depth_spn.visible = false					
			Render_Integ_Emission_Step_spn.visible = false				
			Render_Integ_Glossy_LMax_spn.visible = false
			Render_Integ_Glossy_Ns_spn.visible = false
			Render_Integ_Glossy_DotR_btn.visible = false
			Render_Integ_Glossy_Kd_spn.visible = false
			Render_Integ_Glossy_Ks_spn.visible = false
			Render_Integ_Glossy_Rough_spn.visible = false				
			Render_Integ_IGI_LightPath_spn.visible = false
			Render_Integ_IGI_LightSet_spn.visible = false
			Render_Integ_IGI_MinDist_spn.visible = false
			Render_Integ_IGI_RRThrs_spn.visible = false
			Render_Integ_IGI_MaxDepth_spn.visible = false
			Render_Integ_IGI_GLimit_spn.visible = false
			Render_Integ_IGI_Gather_spn.visible = false		
			Render_Integ_IrCach_MinWht_spn.visible = true
			Render_Integ_IrCach_MinSpac_spn.visible = true
			Render_Integ_IrCach_MaxSpac_spn.visible = true
			Render_Integ_IrCach_MaxAngle_spn.visible = true
			Render_Integ_IrCach_SpecDepth_spn.visible = true
			Render_Integ_IrCach_indDepth_spn.visible = true
			Render_Integ_IrCach_Samples_spn.visible = true				
			Render_Integ_Path_Depth_spn.visible = false				
			Render_Integ_Photon_Caustic_spn.visible = false
			Render_Integ_Photon_Indirect_spn.visible = false
			Render_Integ_Photon_nUsed_spn.visible = false
			Render_Integ_Photon_SpecDepth_spn.visible = false
			Render_Integ_Photon_MaxPhoDepth_spn.visible = false
			Render_Integ_Photon_FinalGath_btn.visible = false
			Render_Integ_Photon_gatherSamp_spn.visible = false
			Render_Integ_Photon_MaxDist_spn.visible = false
			Render_Integ_Photon_GathAng_spn.visible = false
			Render_Integ_Single_Step_spn.visible = false				
			Render_Integ_Whitted_MaDepth_spn.visible = false					
			)
			9: -- Path Tracing
			(
			Render_Integ_AO_Samples_spn.visible = false
			Render_Integ_AO_MaxDist_spn.visible = false				
			Render_Integ_Diffuse_Max_spn.visible = false
			Render_Integ_Diffuse_Ns_spn.visible = false				
			Render_Integ_Dipole_MaxDepth_spn.visible = false	
			Render_Integ_Dipole_MaxError_spn.visible = false	
			Render_Integ_Dipole_MinDist_spn.visible = false	
			Render_Integ_Direct_Depth_spn.visible = false					
			Render_Integ_Emission_Step_spn.visible = false				
			Render_Integ_Glossy_LMax_spn.visible = false
			Render_Integ_Glossy_Ns_spn.visible = false
			Render_Integ_Glossy_DotR_btn.visible = false
			Render_Integ_Glossy_Kd_spn.visible = false
			Render_Integ_Glossy_Ks_spn.visible = false
			Render_Integ_Glossy_Rough_spn.visible = false				
			Render_Integ_IGI_LightPath_spn.visible = false
			Render_Integ_IGI_LightSet_spn.visible = false
			Render_Integ_IGI_MinDist_spn.visible = false
			Render_Integ_IGI_RRThrs_spn.visible = false
			Render_Integ_IGI_MaxDepth_spn.visible = false
			Render_Integ_IGI_GLimit_spn.visible = false
			Render_Integ_IGI_Gather_spn.visible = false		
			Render_Integ_IrCach_MinWht_spn.visible = false
			Render_Integ_IrCach_MinSpac_spn.visible = false
			Render_Integ_IrCach_MaxSpac_spn.visible = false
			Render_Integ_IrCach_MaxAngle_spn.visible = false
			Render_Integ_IrCach_SpecDepth_spn.visible = false
			Render_Integ_IrCach_indDepth_spn.visible =false
			Render_Integ_IrCach_Samples_spn.visible =false				
			Render_Integ_Path_Depth_spn.visible = true				
			Render_Integ_Photon_Caustic_spn.visible = false
			Render_Integ_Photon_Indirect_spn.visible = false
			Render_Integ_Photon_nUsed_spn.visible = false
			Render_Integ_Photon_SpecDepth_spn.visible = false
			Render_Integ_Photon_MaxPhoDepth_spn.visible = false
			Render_Integ_Photon_FinalGath_btn.visible = false
			Render_Integ_Photon_gatherSamp_spn.visible = false
			Render_Integ_Photon_MaxDist_spn.visible = false
			Render_Integ_Photon_GathAng_spn.visible = false
			Render_Integ_Single_Step_spn.visible = false				
			Render_Integ_Whitted_MaDepth_spn.visible = false					
			)
			10: -- Photon Map
			(
			Render_Integ_AO_Samples_spn.visible = false
			Render_Integ_AO_MaxDist_spn.visible = false				
			Render_Integ_Diffuse_Max_spn.visible = false
			Render_Integ_Diffuse_Ns_spn.visible = false				
			Render_Integ_Dipole_MaxDepth_spn.visible = false	
			Render_Integ_Dipole_MaxError_spn.visible = false	
			Render_Integ_Dipole_MinDist_spn.visible = false	
			Render_Integ_Direct_Depth_spn.visible = false					
			Render_Integ_Emission_Step_spn.visible = false				
			Render_Integ_Glossy_LMax_spn.visible = false
			Render_Integ_Glossy_Ns_spn.visible = false
			Render_Integ_Glossy_DotR_btn.visible = false
			Render_Integ_Glossy_Kd_spn.visible = false
			Render_Integ_Glossy_Ks_spn.visible = false
			Render_Integ_Glossy_Rough_spn.visible = false				
			Render_Integ_IGI_LightPath_spn.visible = false
			Render_Integ_IGI_LightSet_spn.visible = false
			Render_Integ_IGI_MinDist_spn.visible = false
			Render_Integ_IGI_RRThrs_spn.visible = false
			Render_Integ_IGI_MaxDepth_spn.visible = false
			Render_Integ_IGI_GLimit_spn.visible = false
			Render_Integ_IGI_Gather_spn.visible = false		
			Render_Integ_IrCach_MinWht_spn.visible = false
			Render_Integ_IrCach_MinSpac_spn.visible = false
			Render_Integ_IrCach_MaxSpac_spn.visible = false
			Render_Integ_IrCach_MaxAngle_spn.visible = false
			Render_Integ_IrCach_SpecDepth_spn.visible = false
			Render_Integ_IrCach_indDepth_spn.visible =false
			Render_Integ_IrCach_Samples_spn.visible =false				
			Render_Integ_Path_Depth_spn.visible = false				
			Render_Integ_Photon_Caustic_spn.visible = true
			Render_Integ_Photon_Indirect_spn.visible = true
			Render_Integ_Photon_nUsed_spn.visible = true
			Render_Integ_Photon_SpecDepth_spn.visible =true
			Render_Integ_Photon_MaxPhoDepth_spn.visible = true
			Render_Integ_Photon_FinalGath_btn.visible = true
			Render_Integ_Photon_gatherSamp_spn.visible = true
			Render_Integ_Photon_MaxDist_spn.visible = true
			Render_Integ_Photon_GathAng_spn.visible = true
			Render_Integ_Single_Step_spn.visible = false				
			Render_Integ_Whitted_MaDepth_spn.visible = false					
			)
			11: -- Single
			(
			Render_Integ_AO_Samples_spn.visible = false
			Render_Integ_AO_MaxDist_spn.visible = false				
			Render_Integ_Diffuse_Max_spn.visible = false
			Render_Integ_Diffuse_Ns_spn.visible = false				
			Render_Integ_Dipole_MaxDepth_spn.visible = false	
			Render_Integ_Dipole_MaxError_spn.visible = false	
			Render_Integ_Dipole_MinDist_spn.visible = false	
			Render_Integ_Direct_Depth_spn.visible = false					
			Render_Integ_Emission_Step_spn.visible = false				
			Render_Integ_Glossy_LMax_spn.visible = false
			Render_Integ_Glossy_Ns_spn.visible = false
			Render_Integ_Glossy_DotR_btn.visible = false
			Render_Integ_Glossy_Kd_spn.visible = false
			Render_Integ_Glossy_Ks_spn.visible = false
			Render_Integ_Glossy_Rough_spn.visible = false				
			Render_Integ_IGI_LightPath_spn.visible = false
			Render_Integ_IGI_LightSet_spn.visible = false
			Render_Integ_IGI_MinDist_spn.visible = false
			Render_Integ_IGI_RRThrs_spn.visible = false
			Render_Integ_IGI_MaxDepth_spn.visible = false
			Render_Integ_IGI_GLimit_spn.visible = false
			Render_Integ_IGI_Gather_spn.visible = false		
			Render_Integ_IrCach_MinWht_spn.visible = false
			Render_Integ_IrCach_MinSpac_spn.visible = false
			Render_Integ_IrCach_MaxSpac_spn.visible = false
			Render_Integ_IrCach_MaxAngle_spn.visible = false
			Render_Integ_IrCach_SpecDepth_spn.visible = false
			Render_Integ_IrCach_indDepth_spn.visible =false
			Render_Integ_IrCach_Samples_spn.visible =false				
			Render_Integ_Path_Depth_spn.visible = false				
			Render_Integ_Photon_Caustic_spn.visible = false
			Render_Integ_Photon_Indirect_spn.visible = false
			Render_Integ_Photon_nUsed_spn.visible = false
			Render_Integ_Photon_SpecDepth_spn.visible =false
			Render_Integ_Photon_MaxPhoDepth_spn.visible = false
			Render_Integ_Photon_FinalGath_btn.visible = false
			Render_Integ_Photon_gatherSamp_spn.visible = false
			Render_Integ_Photon_MaxDist_spn.visible =false
			Render_Integ_Photon_GathAng_spn.visible = false
			Render_Integ_Single_Step_spn.visible = true				
			Render_Integ_Whitted_MaDepth_spn.visible = false					
			)	
			12: -- Use Probes
			(
			Render_Integ_AO_Samples_spn.visible = false
			Render_Integ_AO_MaxDist_spn.visible = false				
			Render_Integ_Diffuse_Max_spn.visible = false
			Render_Integ_Diffuse_Ns_spn.visible = false				
			Render_Integ_Dipole_MaxDepth_spn.visible = false	
			Render_Integ_Dipole_MaxError_spn.visible = false	
			Render_Integ_Dipole_MinDist_spn.visible = false	
			Render_Integ_Direct_Depth_spn.visible = false					
			Render_Integ_Emission_Step_spn.visible = false				
			Render_Integ_Glossy_LMax_spn.visible = false
			Render_Integ_Glossy_Ns_spn.visible = false
			Render_Integ_Glossy_DotR_btn.visible = false
			Render_Integ_Glossy_Kd_spn.visible = false
			Render_Integ_Glossy_Ks_spn.visible = false
			Render_Integ_Glossy_Rough_spn.visible = false				
			Render_Integ_IGI_LightPath_spn.visible = false
			Render_Integ_IGI_LightSet_spn.visible = false
			Render_Integ_IGI_MinDist_spn.visible = false
			Render_Integ_IGI_RRThrs_spn.visible = false
			Render_Integ_IGI_MaxDepth_spn.visible = false
			Render_Integ_IGI_GLimit_spn.visible = false
			Render_Integ_IGI_Gather_spn.visible = false		
			Render_Integ_IrCach_MinWht_spn.visible = false
			Render_Integ_IrCach_MinSpac_spn.visible = false
			Render_Integ_IrCach_MaxSpac_spn.visible = false
			Render_Integ_IrCach_MaxAngle_spn.visible = false
			Render_Integ_IrCach_SpecDepth_spn.visible = false
			Render_Integ_IrCach_indDepth_spn.visible =false
			Render_Integ_IrCach_Samples_spn.visible =false				
			Render_Integ_Path_Depth_spn.visible = false				
			Render_Integ_Photon_Caustic_spn.visible = false
			Render_Integ_Photon_Indirect_spn.visible = false
			Render_Integ_Photon_nUsed_spn.visible = false
			Render_Integ_Photon_SpecDepth_spn.visible =false
			Render_Integ_Photon_MaxPhoDepth_spn.visible = false
			Render_Integ_Photon_FinalGath_btn.visible = false
			Render_Integ_Photon_gatherSamp_spn.visible = false
			Render_Integ_Photon_MaxDist_spn.visible =false
			Render_Integ_Photon_GathAng_spn.visible = false
			Render_Integ_Single_Step_spn.visible = false				
			Render_Integ_Whitted_MaDepth_spn.visible = false					
			)	
			13: -- Whitted
			(
			Render_Integ_AO_Samples_spn.visible = false
			Render_Integ_AO_MaxDist_spn.visible = false				
			Render_Integ_Diffuse_Max_spn.visible = false
			Render_Integ_Diffuse_Ns_spn.visible = false				
			Render_Integ_Dipole_MaxDepth_spn.visible = false	
			Render_Integ_Dipole_MaxError_spn.visible = false	
			Render_Integ_Dipole_MinDist_spn.visible = false	
			Render_Integ_Direct_Depth_spn.visible = false					
			Render_Integ_Emission_Step_spn.visible = false				
			Render_Integ_Glossy_LMax_spn.visible = false
			Render_Integ_Glossy_Ns_spn.visible = false
			Render_Integ_Glossy_DotR_btn.visible = false
			Render_Integ_Glossy_Kd_spn.visible = false
			Render_Integ_Glossy_Ks_spn.visible = false
			Render_Integ_Glossy_Rough_spn.visible = false				
			Render_Integ_IGI_LightPath_spn.visible = false
			Render_Integ_IGI_LightSet_spn.visible = false
			Render_Integ_IGI_MinDist_spn.visible = false
			Render_Integ_IGI_RRThrs_spn.visible = false
			Render_Integ_IGI_MaxDepth_spn.visible = false
			Render_Integ_IGI_GLimit_spn.visible = false
			Render_Integ_IGI_Gather_spn.visible = false		
			Render_Integ_IrCach_MinWht_spn.visible = false
			Render_Integ_IrCach_MinSpac_spn.visible = false
			Render_Integ_IrCach_MaxSpac_spn.visible = false
			Render_Integ_IrCach_MaxAngle_spn.visible = false
			Render_Integ_IrCach_SpecDepth_spn.visible = false
			Render_Integ_IrCach_indDepth_spn.visible =false
			Render_Integ_IrCach_Samples_spn.visible =false				
			Render_Integ_Path_Depth_spn.visible = false				
			Render_Integ_Photon_Caustic_spn.visible = false
			Render_Integ_Photon_Indirect_spn.visible = false
			Render_Integ_Photon_nUsed_spn.visible = false
			Render_Integ_Photon_SpecDepth_spn.visible =false
			Render_Integ_Photon_MaxPhoDepth_spn.visible = false
			Render_Integ_Photon_FinalGath_btn.visible = false
			Render_Integ_Photon_gatherSamp_spn.visible = false
			Render_Integ_Photon_MaxDist_spn.visible =false
			Render_Integ_Photon_GathAng_spn.visible = false
			Render_Integ_Single_Step_spn.visible = false				
			Render_Integ_Whitted_MaDepth_spn.visible = true					
			)				
		)
	)
	
	on Render_Filter_ddl selected arg do (
		case Render_Filter_ddl.selection of (
			1: -- Box
			(
			Render_Filter_Box_X.visible = true
			Render_Filter_Box_Y.visible = true
			Render_Filter_Gauss_X.visible = false
			Render_Filter_Gauss_Y.visible = false
			Render_Filter_Gauss_Alpah.visible = false
			Render_Filter_Mitch_X.visible = false
			Render_Filter_Mitch_Y.visible = false
			Render_Filter_Mitch_Sharp.visible = false
			Render_Filter_Sinc_X.visible = false
			Render_Filter_Sinc_Y.visible = false
			Render_Filter_Sinc_Tau.visible = false
			Render_Filter_Tri_X.visible = false
			Render_Filter_Tri_Y.visible = false				
			)
			2: -- Gaussian
			(
			Render_Filter_Box_X.visible = false
			Render_Filter_Box_Y.visible = false
			Render_Filter_Gauss_X.visible = true
			Render_Filter_Gauss_Y.visible = true
			Render_Filter_Gauss_Alpah.visible = true
			Render_Filter_Mitch_X.visible = false
			Render_Filter_Mitch_Y.visible = false
			Render_Filter_Mitch_Sharp.visible = false
			Render_Filter_Sinc_X.visible = false
			Render_Filter_Sinc_Y.visible = false
			Render_Filter_Sinc_Tau.visible = false
			Render_Filter_Tri_X.visible = false
			Render_Filter_Tri_Y.visible = false					
			)
			3: -- Mitchell
			(
			Render_Filter_Box_X.visible = false
			Render_Filter_Box_Y.visible = false
			Render_Filter_Gauss_X.visible = false
			Render_Filter_Gauss_Y.visible = false
			Render_Filter_Gauss_Alpah.visible = false
			Render_Filter_Mitch_X.visible =true
			Render_Filter_Mitch_Y.visible = true
			Render_Filter_Mitch_Sharp.visible = true
			Render_Filter_Sinc_X.visible = false
			Render_Filter_Sinc_Y.visible = false
			Render_Filter_Sinc_Tau.visible = false
			Render_Filter_Tri_X.visible = false
			Render_Filter_Tri_Y.visible = false				
			)
			4: -- Sinc
			(
			Render_Filter_Box_X.visible = false
			Render_Filter_Box_Y.visible = false
			Render_Filter_Gauss_X.visible = false
			Render_Filter_Gauss_Y.visible = false
			Render_Filter_Gauss_Alpah.visible = false
			Render_Filter_Mitch_X.visible =false
			Render_Filter_Mitch_Y.visible = false
			Render_Filter_Mitch_Sharp.visible = false
			Render_Filter_Sinc_X.visible = true
			Render_Filter_Sinc_Y.visible = true
			Render_Filter_Sinc_Tau.visible = true
			Render_Filter_Tri_X.visible = false
			Render_Filter_Tri_Y.visible = false					
			)
			5: -- Triangle
			(
			Render_Filter_Box_X.visible = false
			Render_Filter_Box_Y.visible = false
			Render_Filter_Gauss_X.visible = false
			Render_Filter_Gauss_Y.visible = false
			Render_Filter_Gauss_Alpah.visible = false
			Render_Filter_Mitch_X.visible =false
			Render_Filter_Mitch_Y.visible = false
			Render_Filter_Mitch_Sharp.visible = false
			Render_Filter_Sinc_X.visible = false
			Render_Filter_Sinc_Y.visible = false
			Render_Filter_Sinc_Tau.visible =false
			Render_Filter_Tri_X.visible = true
			Render_Filter_Tri_Y.visible = true				
			)
		)
	)
	)

-- Camera
rollout PBRT_Camera "Camera" 
	(
	dropDownList Render_Camera_ddl "" pos:[14,12] width:122 height:21 items:#("Enviorment", "Orthographic", "Camera", "Perspective") selection:3
	
	spinner 		PBRT_Camera_ShOpen_spn 	"Shutter Open:" 	pos:[55,42] width:82 height:16 range:[0,10000,0]
	spinner 		PBRT_Camera_ShClose_spn 	"Shutter Close:" 	pos:[55,64] width:82 height:16 range:[0,10000,1]
	spinner 		PBRT_Camera_LensR_spn 		"Lens Radius:" 	pos:[55,87] width:82 height:16 range:[0,10000,0]
	spinner 		PBRT_Camera_FDist_spn 		"Focal Distance:" pos:[55,110] width:82 height:16 range:[0,1e+007,1e+006]
	listBox 			PBRT_Camera_Cameras_lbx 	"Cameras:" 			pos:[150,12] width:150 height:9

	Local theCameras = #()
	
	on PBRT_Camera open do (
		for o in cameras where classof o != Targetobject do  append theCameras o.name
		PBRT_Camera_Cameras_lbx.items = theCameras
		if theCameras.count != 0 do CurrentCamera = (getNodeByName theCameras[1])
		
		--update the values from the camera
		if theCam == undefined then (
			Render_Camera_ddl.selection = 2
		)else (
				-- we need to change it to  trget camera to calculate the distance

			theCam.type = #target
			theDist = distance theCam.pos theCam.target.pos
			theCam.type = #free
			PBRT_Camera_FDist_spn.value = theDist		
			
			theCam = getNodeByName PBRT_Camera_Cameras_lbx.items[1]			
		)

	)
	
	on PBRT_Camera_Cameras_lbx selected index do (
		CurrentCamera =  (getNodeByName PBRT_Camera_Cameras_lbx.items[index])
	)
	
	
	)

-- Tools
rollout PBRT_Tools "Tools" 
	(
		checkbox 	uncPath_chk	"Convert to UNC Paths"  	width:130 height:17 offset:[0,10] 

		edittext	uncPath_etxt	"Path:"  							width:260 height:17 across:2 --offset:[0,10]	
		button 		uncPath_btn 	"..." 		 						width:15 height:17 align:#right
		
		button		CopyBitmaps	"Copy Maps"	width:70 height:18 align:#left tooltip:"Click to select destination folder and copy all maps to that folder."
		
		GroupBox 	Tools_Bitmap_Grp "Bitmaps" 	pos:[5,0] 		width:293 	height:85 

		on PBRT_Tools open do (
			uncPath_etxt.text = "textures"
		)
		
		Local theMapsPath = undefined 
		on uncPath_btn pressed do (

			SAVEPATH = getSavePath caption:"-- Gtools --  Select Path where you want the files to be copied." --initialDir:"d:"
				if (SAVEPATH  != undefined) then
				(			
					theMapsPath = SAVEPATH
					uncPath_etxt.text = SAVEPATH --"\.. " + (OrifilePathReduced as string)	
				)
			)
			
		on CopyBitmaps pressed do (
			CopyMaps uncPath_etxt.text
		)
	)
	------------
--closeRolloutFloater PBRTMainFloater
PBRTMainFloater = newRolloutFloater "Gtools PBRT Render Exporter v0.50" 315 780

addRollout PBRT_System PBRTMainFloater rolledUp:true
addRollout PBRT_Output PBRTMainFloater
addRollout PBRT_Render PBRTMainFloater
addRollout PBRT_Camera PBRTMainFloater
addRollout PBRT_Tools PBRTMainFloater rolledUp:true



)
