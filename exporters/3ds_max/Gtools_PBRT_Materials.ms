/*


Known Issues:
- opacity its not shown on the previews, as pbrt dosnt support opacity on the primitives, and i'm using a sphere
*/
-- PBRT Materials
Global PBRT_EXE_PATH = "C:\\Program Files\\Graphics\\Autodesk\\3ds Max 2010\\export\\"
	

	fn PBRT_clampRange 
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
		
	-- function taken from Diffuse Mapper  by André Hotz , http://www.scriptspot.com/andrehotz/
	fn PBRT_FindAveColor Dif =
		(
			difmap = Dif
			cbitmap = rendermap difmap size:[300, 300]

			bh = cbitmap.height
			bw = cbitmap.width
			numpix = bh * bw
			
			redline   = 0.0
			greenline = 0.0
			blueline  = 0.0
			
			for line = 1 to bh do
			(
				colorRow = getPixels cbitmap [0,(line-1)] bw
				for c=1 to bw do
				(
					redLine   = redLine   + colorRow[c].red
					greenLine = greenLine + colorRow[c].green
					blueLine  = blueLine  + colorRow[c].blue
				)
				colorRow = #()
			)
			avgRed   = RedLine   / numPix
			avgGreen = GreenLine / numPix
			avgBlue  = BlueLine  / numPix
			avgColor = color avgRed avgGreen avgBlue

			return avgColor
			close cbitmap
		)	
	
	fn PBRT_getMaterialList gMatList =
			(
				matList = #()
				matClass = classof gmatList
				if matClass == MultiMaterial then
				(
					for sm = 1 to gmatList.count do
					(
					append matList gmatList[sm]
					)
				)else
				append matList gmatList
				
				return matList 
			)
	
	fn PbrtColor c =
		(
			return ( ((c.r / 255)) as string + " " + ((c.g / 255)) as string + " " + ((c.b / 255)) as string )
		)
	
	fn PbrtAddBack str =
		(
			if (str == undefined) do 
				str = ""
			
			str = substituteString str "\\" "\\\\"
			
			return ("\"" + str + "\"") 
		)
	
	fn Pbrt_Textures type &out_file mat =
		(		
			Local map --= obj.material.KdTexture
			Local ExportedTex = #()
			
			Local textureString = #("","")

			Local MapType

			if type ==	"KdTexture" then (
				MapType = "imagemap" 
				map = mat.KdTexture
			)
				else if type == "bumpMap" then (
					MapType = "scale"
					map = mat.bumpMap
				) else if type == "sigma" then (
						MapType = "sigma"
						map = mat.sigmaMap				
					)else if type == "KsMap" then (
						MapType = "KsMap"
						map = mat.KsTexture				
						)else if type == "RoughMap" then (
							MapType = "RoughMap"
							map = mat.RoughnessTexture				
						)else if type == "KrMap" then (
								MapType = "KrMap"
								map = mat.KrTexture				
							)else if type == "OpacityMap" then (
								MapType = "OpacityMap"
								map = mat.opacityTexture				
							)else if type == "LightMap" then (
								MapType = "LightMap"
								map = mat.LightCMap				
							)else if type == "AlphaMap" then (
								MapType = "AlphaMap"
								map = mat.alphaImage				
							)

					
			if (classof map == PBRT_BitMap) then
				(
		
					local wrap 			= #("repeat", "black", "clamp") 	
					local mapType	= #("uv", "spherical", "cylindrical", "planar")
					local mapTypeID	= map.uv
					local wrapID 		= map.wrap
					local tri				= map.tri
					local anist			= map.anist
					local gamma 		= map.gamma	
					local file			= PbrtAddBack map.PBRT_BitMap.filename

					if type == "bumpMap" or type == "sigma"  or type == "AlphaMap" then (
						format "Texture \"%\" \"float\" \"imagemap\"\n" map.name to:out_file
					)else (
					format "Texture \"%\" \"color\" \"imagemap\"\n" map.name to:out_file
					)
					format "\t \"string mapping\" \"%\"\n" mapType[mapTypeID] to:out_file			
					format "\t \"string wrap\" [\"%\"] \n" wrap[wrapID] to:out_file
					format "\t \"string filename\" [%]\n"  file to:out_file
					format "\t \"float gamma\" [%]\n" gamma to:out_file
					format "\t \"string mapping\" [\"%\"]\n" mapType[mapTypeID] to:out_file
					format "\t \"float uscale\" [%]\n" (map.uscale * 2) to:out_file
					format "\t \"float vscale\" [%]\n" (map.vscale * 2) to:out_file
					format "\t \"float udelta\" [%]\n" map.udelta to:out_file
					format "\t \"float vdelta\" [%]\n" map.vdelta to:out_file	

					if type == "bumpMap" then (
						Local MapName = map.name + "_bump"
						format "Texture \"%\" \"float\" \"scale\" \"texture tex1\" \"%\" \"float tex2\" [%]  \n" MapName map.name mat.bumpMapAmount to:out_file
					)
					
					if type == "sigma" then (
						Local MapName = map.name + "_sigma"
						format "Texture \"%\" \"float\" \"scale\" \"texture tex1\" [\"%\"] \"float tex2\" [%]  \n" MapName map.name mat.sigma to:out_file
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
						format "Texture \"%\" \"color\" \"scale\" \"texture tex1\" [\"%\"] \"float tex2\" [% % %]  \n" MapName map.name 1.0 1.0 1.0 to:out_file
					)
					
					if type == "RoughMap" then (
						Local MapName = map.name + "_RoughMap"
						format "Texture \"%\" \"color\" \"scale\" \"texture tex1\" [\"%\"] \"float tex2\" [%]  \n" MapName map.name mat.roughness to:out_file
					)
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
						format "Texture \"%\" \"scale\" \"texture tex1\" \"%\" \"float tex2\" [%]  \n" MapName map.name mat.bumpMapAmount to:out_file
					)			
					
	--				format "Material \"matte\" \"texture Kd\" [\"%\"]\n"  map.name to:out_file			
			)
			else if (classof map == PBRT_ConstantMap) then (

					format "Texture \"%\" \"color\" \"scale\"\n" map.name  to:out_file
					format "\t \"color tex1\" [%] \n" (PbrtColor map.color)  to:out_file
					format "\t \"color tex2\" [%] \n" (PbrtColor this.kd)  to:out_file			
					
	--				format "Material \"matte\" \"texture Kd\" [\"%\"]\n"  map.name to:out_file			
			)		
		)	
	
	fn Pbrt_Material &out_file mat = (

		-- i'm not sure why i cant test for the material class???  so i had to make a string from it.
		theClass = (classof mat) as string

			-- write object name
			format "\n### Object: % ### \n\n" mat.name  to:out_file			

			-- ************     Matte Material     ***********************
							if theClass == "PBRT_Matte" then (
						
								local Kd = PbrtColor mat.Kd
				
								Local hasDiffuseMap = false
								Local hasSigmaMap	= false					
								Local hasBumpMap	= false
								Local hasAlphaMap	= false
	--							Local hasLightMap	= false
								
								------- *** Create the Textures ***
									
										-- Diffuse Texture
										if mat.KdTextureEnable == true and mat.KdTexture != undefined then (
											hasDiffuseMap = true

											Pbrt_Textures "KdTexture" out_file mat	-- create the texture
										)
										-- Sigma Texture
										if mat.sigmaMapEnabled == true and mat.sigmaMap != undefined then (
											hasSigmaMap = true

												Pbrt_Textures "sigma" out_file mat -- create the texture
										)					
										-- Bump Map Texture
										if mat.bumpMapEnable == true and mat.bumpMap != undefined then (
											hasBumpMap	= true
											
													Pbrt_Textures "bumpMap" out_file mat -- create the texture		
										)
										-- Alpha Map Texture
										if mat.alphaEnable == true and mat.alphaImage != undefined then (
											hasAlphaMap	= true
											
													Pbrt_Textures "AlphaMap"	out_file mat -- create the texture		
										)																		
																	
								-- *** Create The Material *** --

								format "MakeNamedMaterial \"%\" \n" mat.name  to:out_file
								format "\t \"string type\" [\"matte\"] \n" to:out_file		
								
								
								-- if it dosent ha a map or its not beenuse then use the color
								if mat.KdTextureEnable == false  or mat.KdTexture == undefined  then (
								format "\t \"color Kd\" [%]\n"  kd to:out_file
								)
															
								if mat.sigma > 0.0 then format "\t \"float sigma\" [%]\n"  mat.sigma to:out_file
							
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
							
								------- Create the Textures
								if mat.KdTextureEnable == true and mat.KdTexture != undefined then (
									hasDiffuseMap = true

										Pbrt_Textures "KdTexture" out_file	mat -- create the texture
									)
								
								if mat.KsTextureEnable == true and mat.KsTexture != undefined then (
									hasKsMap = true
									
											Pbrt_Textures "KsMap" out_file	mat -- create the texture
									)					
								
								if mat.bumpMapEnable == true and mat.bumpMap != undefined then (
									hasBumpMap	= true

											Pbrt_Textures "bumpMap"	out_file mat -- create the texture
									)
									
								if mat.RoughnessTextureEnable == true and mat.RoughnessTexture != undefined then (
									hasRoughMap	= true
									
												Pbrt_Textures "RoughMap" out_file	mat -- create the texture
									)	
								
								-- Alpha Map Texture
								if mat.alphaEnable == true and mat.alphaImage != undefined then ( --**AQUI
										hasAlphaMap	= true

												Pbrt_Textures "AlphaMap"	out_file mat -- create the texture
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
						if mat.KdTextureEnable == true and mat.KdTexture != undefined then (
							hasDiffuseMap = true
							
							Pbrt_Textures "KdTexture" out_file	mat -- create the texture
						)
						
						if mat.KsTextureEnable == true and mat.KsTexture != undefined then (
							hasKsMap = true

							Pbrt_Textures "KsMap" out_file mat	-- create the texture
						)					
						
						if mat.KrTextureEnable == true and mat.KrTexture != undefined then (
							hasKrMap	= true

							Pbrt_Textures "KrMap" out_file mat	-- create the texture
						)							
						
						if mat.opacityTextureEnable == true and mat.opacityTexture != undefined then (
							hasOpacityMap	= true

							Pbrt_Textures "OpacityMap" out_file mat -- create the texture
						)
						
						if mat.bumpMapEnable == true and mat.bumpMap != undefined then (
							hasBumpMap	= true
																																		
							Pbrt_Textures "bumpMap" out_file mat	-- create the texture
						)
							
						if mat.RoughnessTextureEnable == true and mat.RoughnessTexture != undefined then (
							hasRoughMap	= true
							
							Pbrt_Textures "RoughMap" out_file	mat -- create the texture
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
					if mat.bumpMapEnable == true and mat.bumpMap != undefined then (
								hasBumpMap	= true

								Pbrt_Textures  "bumpMap" out_file	mat -- create the texture
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
					

						------- Create the Textures
						if mat.KdTextureEnable == true and mat.KdTexture != undefined then (
							hasDiffuseMap = true

								Pbrt_Textures "KdTexture" out_file	mat -- create the texture
							)
						
						if mat.KsTextureEnable == true and mat.KsTexture != undefined then (
							hasKsMap = true

								Pbrt_Textures "KsMap" out_file mat -- create the texture
							)					
						
						if mat.bumpMapEnable == true and mat.bumpMap != undefined then (
							hasBumpMap	= true

								Pbrt_Textures "bumpMap"	out_file mat -- create the texture
							)
							
						if mat.RoughnessTextureEnable == true and mat.RoughnessTexture != undefined then (
							hasRoughMap	= true

								Pbrt_Textures "RoughMap" out_file	mat -- create the texture
							)	
							
						if mat.reflectTextureEnable == true and mat.reflectTexture != undefined then (
							hasReflectMap	= true

								Pbrt_Textures "ReflectMap" out_file	mat -- create the texture
							)							
							
						if mat.transmitTextureEnable == true and mat.transmitTexture != undefined then (
							hasTransmitMap	= true

								Pbrt_Textures "TransmitMap" out_file	mat -- create the texture
							)	
						
						-- Alpha Map Texture
						if mat.alphaEnable == true and mat.alphaImage != undefined then ( --**AQUI
								hasAlphaMap	= true

									Pbrt_Textures "AlphaMap"	out_file mat -- create the texture
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
					
						------- Create the Textures
						if mat.KdTextureEnable == true and mat.KdTexture != undefined then (
							hasDiffuseMap = true

								Pbrt_Textures "KdTexture" out_file	mat -- create the texture
							)
						
						if mat.KsTextureEnable == true and mat.KsTexture != undefined then (
							hasKsMap = true

								Pbrt_Textures "KsMap" out_file mat 	-- create the texture
							)					
						
						if mat.bumpMapEnable == true and mat.bumpMap != undefined then (
							hasBumpMap	= true

								Pbrt_Textures "bumpMap"	out_file mat -- create the texture
							)
							
						if mat.uRoughnessTextureEnable == true and mat.uRoughnessTexture != undefined then (
							hasURoughMap	= true

								Pbrt_Textures "uRoughMap" out_file	 mat -- create the texture
							)	
							
						if mat.vRoughnessTextureEnable == true and mat.vRoughnessTexture != undefined then (
							hasVRoughMap	= true

								Pbrt_Textures "vRoughMap" out_file	 mat -- create the texture
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
					

						------- Create the Textures				
						
						if mat.bumpMapEnable == true and mat.bumpMap != undefined then (
							hasBumpMap	= true
							
								Pbrt_Textures "bumpMap"	out_file mat -- create the texture
							)
							
						if mat.RoughnessTextureEnable == true and mat.RoughnessTexture != undefined then (
							hasRoughMap	= true

								Pbrt_Textures "RoughMap" out_file	mat -- create the texture
							)										

					
					-- *** Create The Material *** --
					format "MakeNamedMaterial \"%\" \n" mat.name  to:out_file
					format "\t \"string type\" [\"metal\"] \n" to:out_file
					
					if mat.RoughnessTextureEnable == false then format "\t \"float roughness\"  [%] \n"  mat.Roughness to:out_file

										
						format "\t \"spectrum eta\" % \n" (PbrtAddBack mat.eta) to:out_file
						format "\t \"spectrum k\" % \n" (PbrtAddBack mat.k) to:out_file

						if hasRoughMap then (
							Local RoughMapName = mat.RoughnessTexture.name + "_RoughMap"						
							format "\t \"texture roughness\" \"%\"\n"  RoughMapName to:out_file								
						)
						if hasBumpMap then (
							Local BumpMapName = mat.bumpMap.name + "_bump"						
							format "\t \"texture bumpmap\" [\"%\"]\n"  BumpMapName to:out_file							
						)				
									
				)
				
							format "\n### End Object: % ### \n" mat.name  to:out_file				
	)

	fn PBRT_MakePreview mat RendMode PrevQuality = (
			
		Local Res = 140
		Local fname = "preview"
		Local SAVEPATH =  (GetDir #export) + "\\" --+ fname

	
		Local mainTexFile = SAVEPATH + fname + ".pbrt"

		Local out_file = createfile mainTexFile


		-- Create the Scene
		
		format "# pbrt v2.0 Scene File Ver % \n" ver to:out_file
		format "# Exported by GtoolExporter \n\n" to:out_file
		
--		format "LookAt  0.18423 -0.256156 0.238845  -0.305277 0.327216 -0.409275  -0.416603 0.496489 0.761538 \n"to:out_file
		format "LookAt  0.205405 -0.288824 0.173344  -0.383984 0.413583 -0.225712  -0.256508 0.305694 0.916927 \n"to:out_file
		format "Camera \"perspective\" \"float fov\" [43.0236]  \"float shutteropen\" [0]   \"float shutterclose\" [1] \n"  to:out_file
		format "Film \"image\" \"integer xresolution\"  [%]  \"integer yresolution\"  [%]   \"string filename\"  [\"preview.exr\"]\n"  Res Res  to:out_file
		format "PixelFilter \"mitchell\" \n"  to:out_file

		
		if RendMode == 1 do ( -- Direct Lighting	
			if PrevQuality == 1 do (
				format "Sampler \"lowdiscrepancy\" \"integer pixelsamples\" [4] \n"  to:out_file
				format "SurfaceIntegrator \"directlighting\" \"integer maxdepth\" [3] "  to:out_file
			)
			if PrevQuality == 2 do (
				format "Sampler \"lowdiscrepancy\" \"integer pixelsamples\" [8] \n"  to:out_file
				format "SurfaceIntegrator \"directlighting\" \"integer maxdepth\" [6] "  to:out_file			
			)
			if PrevQuality == 3 do (
				format "Sampler \"lowdiscrepancy\" \"integer pixelsamples\" [16] \n"  to:out_file
				format "SurfaceIntegrator \"directlighting\" \"integer maxdepth\" [8] "  to:out_file				
			)	
			if PrevQuality == 4 do (
				format "Sampler \"lowdiscrepancy\" \"integer pixelsamples\" [32] \n"  to:out_file
				format "SurfaceIntegrator \"directlighting\" \"integer maxdepth\" [8] "  to:out_file				
			)		
			format "Renderer \"sample\" \"float stepsize\" [1.0]  \n"  to:out_file
		)
		
		if RendMode == 2 do ( -- Metropolis
			format "SurfaceIntegrator \"path\" \"integer maxdepth\" [3] \n"  to:out_file			
			
			if PrevQuality == 1 do (
				format "Renderer \"metropolis\"  \"float largestepprobability\" [0.5]  \"integer samplesperpixel\" [5]  \"integer bootstrapsamples\" [10]"  to:out_file
				format "\"integer directsamples\" [3]  \"integer maxconsecutiverejects\" [10] \"integer maxdepth\" [3] \n"  to:out_file
			)				
			if PrevQuality == 2 do (
				format "Renderer \"metropolis\"  \"float largestepprobability\" [0.5]  \"integer samplesperpixel\" [32]  \"integer bootstrapsamples\" [10]"  to:out_file
				format "\"integer directsamples\" [6]  \"integer maxconsecutiverejects\" [10] \"integer maxdepth\" [4] \n"  to:out_file
			)			
			if PrevQuality == 3 do (
				format "Renderer \"metropolis\"  \"float largestepprobability\" [0.5]  \"integer samplesperpixel\" [64]  \"integer bootstrapsamples\" [10]"  to:out_file
				format "\"integer directsamples\" [14]  \"integer maxconsecutiverejects\" [10] \"integer maxdepth\" [6] \n"  to:out_file
			)			
			if PrevQuality == 4 do (
				format "Renderer \"metropolis\"  \"float largestepprobability\" [0.5]  \"integer samplesperpixel\" [248]  \"integer bootstrapsamples\" [10]"  to:out_file
				format "\"integer directsamples\" [32]  \"integer maxconsecutiverejects\" [10] \"integer maxdepth\" [6] \n"  to:out_file
			)
		)
		
		format "Accelerator \"bvh\"\n"  to:out_file


		format "WorldBegin\n"  to:out_file
		-- Enviorment 
		format "AttributeBegin\n"  to:out_file
		format "\t LightSource \"infinite\"  \"color L\" [0.9 0.9 0.9] \"color scale\" [1.0 1.0 1.0] \n"  to:out_file		
		if PrevQuality == 1 do format "\"integer nsamples\" [4] " to:out_file						
		if PrevQuality == 2 do format "\"integer nsamples\" [8] " to:out_file
		if PrevQuality == 3 do format "\"integer nsamples\" [16] " to:out_file
		if PrevQuality == 4 do format "\"integer nsamples\" [32] " to:out_file			
		format "AttributeEnd \n"  to:out_file
/*	*/
		-- ***** Lights ***** 
		format "TransformBegin \n"  to:out_file
		format "\t Translate -0.5 -0.6 0.55 \n"  to:out_file
--		format "\t Translate -0.6 -0.801089 0.6 \n"  to:out_file				
		format "\t LightSource \"spot\"  \"point from\" [0 0 0]  \"point to\" [0.577548 0.409736 -0.6] "  to:out_file		
--		format "\t LightSource \"spot\"  \"point from\" [0 0 0]  \"point to\" [0.0719 -0.1244 0.0061] "  to:out_file
		format "\"color I\" [.9 .9 .9]  \"float coneangle\" [27.4] \"float conedeltaangle\" [17.6]  \n"  to:out_file
		format "TransformEnd \n"  to:out_file

		
		-- Materials
		-- Floor 
		format "Texture \"CheckerMap\" \"color\" \"checkerboard\" \"string mapping\" \"uv\"  \"color tex1\" [0.7 0.7 0.7]"  to:out_file
		format "\"color tex2\" [0.85 0.85 0.85]   \"string aamode\" [\"closedform\"] \"float uscale\" [8.0] \"float vscale\" [8.0]\n"  to:out_file
		format "MakeNamedMaterial \"Piso\" \"string type\" [\"matte\"]  \"texture Kd\" [\"CheckerMap\"]\n"  to:out_file
		-- Sphere

		Pbrt_Material out_file mat
		-- Geometry
/*			
		-- Floor
		format "AttributeBegin \n" to:out_file
		format "NamedMaterial \"Piso\"  \n" to:out_file
		format "\t Shape \"trianglemesh\" \"integer indices\" \n" to:out_file
		format "\t [0 1 2  3 4 5] \n" to:out_file
		format "\t \"point P\"  \n" to:out_file
		format "\t [-5.0 5.0 0.0  -5.0 -5.0 0.0  5.0 5.0 0.0  5.0 -5.0 0.0  5.0 5.0 0.0  -5.0 -5.0 0.0] \n" to:out_file		
		format "\t \"float uv\"  \n" to:out_file	
		format "\t [-4.495 5.495  -4.495 -4.495  5.495 5.495  5.495 -4.495  5.495 5.495  -4.495 -4.495] \n" to:out_file
		format "\t \"normal N\"  \n" to:out_file	
		format "\t [0.0 0.0 1.0  0.0 0.0 1.0  0.0 0.0 1.0  0.0 0.0 1.0  0.0 0.0 1.0  0.0 0.0 1.0] \n" to:out_file
		format "AttributeEnd \n" to:out_file	

		-- Sphere
	
		format "AttributeBegin \n" to:out_file
		format "NamedMaterial \"%\"  \n" mat.name to:out_file
		format "\t Transform [1.0 0.0 0.0 0  0.0 1.0 0.0 0   0.0 0.0 1.0 0    0.079721 -0.12887 0.100483 1]     \n" to:out_file
		format "\t Shape \"sphere\" \"float radius\" [0.1] \n" to:out_file	
		format "AttributeEnd \n" to:out_file	
*/

		-- **** Object: Piso **** 
		format "AttributeBegin \n" to:out_file 
		format "NamedMaterial \"Piso\"\n" to:out_file 

		format "Shape \"trianglemesh\" \"integer indices\"\n"  to:out_file 
		format "[0 1 2	3 4 5]\n"  to:out_file 
		format "\"point P\" [-5.0 5.0 0.0	-5.0 -5.0 0.0	5.0 5.0 0.0	5.0 -5.0 0.0	5.0 5.0 0.0	-5.0 -5.0 0.0]\n"  to:out_file 
		format "\"float uv\" [-5.96317 6.96317	-5.96317 -5.96317	6.96317 6.96317	6.96317 -5.96317	6.96317 6.96317	-5.96317 -5.96317]\n"  to:out_file 
		format "\"normal N\" [0.0 0.0 1.0	0.0 0.0 1.0	0.0 0.0 1.0	0.0 0.0 1.0	0.0 0.0 1.0	0.0 0.0 1.0]\n"  to:out_file 
		format "AttributeEnd \n" to:out_file
		--**** End Object: Piso **** 
		

		--**** Object: Sphere **** 
		format "AttributeBegin \n" to:out_file
		format "NamedMaterial \"%\"  \n" mat.name to:out_file
		
		-- Geometry Area Lights
		if classof mat == PBRT_Light or \
			(classof mat== PBRT_Matte and mat.LightColorEnable) or \
				(classof mat == PBRT_Uber and mat.LightColorEnable) then (
					

					
			theColor = (PbrtColor mat.LightColor)
			theColor2 = (PbrtColor (mat.LightColor * 5))

			format "\t AreaLightSource \"area\" \"color L\" [%]\n" theColor to:out_file
			format "\t \"color scale\" [% % %]\n" mat.lightMult mat.lightMult mat.lightMult to:out_file	
			format "\t \"integer nsamples\" [%] \n" mat.LightSamples to:out_file
		
			Local theName = uniquename (mat.name + "_Mat")
			
			format "\t Texture \"%\" \"color\" \"scale\"\n" theName to:out_file
			format "\t \"color tex1\" [%] \n" theColor2 to:out_file
--					format "\t \"color tex2\" [%] \n" theColor  to:out_file	
			format "\n" 	to:out_file						
			
		)
				
		format "Shape \"trianglemesh\" \"integer indices\" [0 1 2	3 4 5	6 7 8	9 10 11	12 13 14	15 16 17	18 19 20	21 22 23	24 25 26	27 28 29	30 31 32	33 34 35	36 37 38	39 40 41	42 43 44	45 46 47	48 49 50	51 52 53	54 55 56	57 58 59	60 61 62	63 64 65	66 67 68	69 70 71	72 73 74	75 76 77	78 79 80	81 82 83	84 85 86	87 88 89	90 91 92	93 94 95	96 97 98	99 100 101	102 103 104	105 106 107	108 109 110	111 112 113	114 115 116	117 118 119	120 121 122	123 124 125	126 127 128	129 130 131	132 133 134	135 136 137	138 139 140	141 142 143	144 145 146	147 148 149	150 151 152	153 154 155	156 157 158	159 160 161	162 163 164	165 166 167	168 169 170	171 172 173	174 175 176	177 178 179	180 181 182	183 184 185	186 187 188	189 190 191	192 193 194	195 196 197	198 199 200	201 202 203	204 205 206	207 208 209	210 211 212	213 214 215	216 217 218	219 220 221	222 223 224	225 226 227	228 229 230	231 232 233	234 235 236	237 238 239	240 241 242	243 244 245	246 247 248	249 250 251	252 253 254	255 256 257	258 259 260	261 262 263	264 265 266	267 268 269	270 271 272	273 274 275	276 277 278	279 280 281	282 283 284	285 286 287	288 289 290	291 292 293	294 295 296	297 298 299	300 301 302	303 304 305	306 307 308	309 310 311	312 313 314	315 316 317	318 319 320	321 322 323	324 325 326	327 328 329	330 331 332	333 334 335	336 337 338	339 340 341	342 343 344	345 346 347	348 349 350	351 352 353	354 355 356	357 358 359	360 361 362	363 364 365	366 367 368	369 370 371	372 373 374	375 376 377	378 379 380	381 382 383	384 385 386	387 388 389	390 391 392	393 394 395	396 397 398	399 400 401	402 403 404	405 406 407	408 409 410	411 412 413	414 415 416	417 418 419	420 421 422	423 424 425	426 427 428	429 430 431	432 433 434	435 436 437	438 439 440	441 442 443	444 445 446	447 448 449	450 451 452	453 454 455	456 457 458	459 460 461	462 463 464	465 466 467	468 469 470	471 472 473	474 475 476	477 478 479	480 481 482	483 484 485	486 487 488	489 490 491	492 493 494	495 496 497	498 499 500	501 502 503	504 505 506	507 508 509	510 511 512	513 514 515	516 517 518	519 520 521	522 523 524	525 526 527	528 529 530	531 532 533	534 535 536	537 538 539	540 541 542	543 544 545	546 547 548	549 550 551	552 553 554	555 556 557	558 559 560	561 562 563	564 565 566	567 568 569	570 571 572	573 574 575	576 577 578	579 580 581	582 583 584	585 586 587	588 589 590	591 592 593	594 595 596	597 598 599	600 601 602	603 604 605	606 607 608	609 610 611	612 613 614	615 616 617	618 619 620	621 622 623	624 625 626	627 628 629	630 631 632	633 634 635	636 637 638	639 640 641	642 643 644	645 646 647	648 649 650	651 652 653	654 655 656	657 658 659	660 661 662	663 664 665	666 667 668	669 670 671	672 673 674	675 676 677	678 679 680	681 682 683	684 685 686	687 688 689	690 691 692	693 694 695	696 697 698	699 700 701	702 703 704	705 706 707	708 709 710	711 712 713	714 715 716	717 718 719	720 721 722	723 724 725	726 727 728	729 730 731	732 733 734	735 736 737	738 739 740	741 742 743	744 745 746	747 748 749	750 751 752	753 754 755	756 757 758	759 760 761	762 763 764	765 766 767	768 769 770	771 772 773	774 775 776	777 778 779	780 781 782	783 784 785	786 787 788	789 790 791	792 793 794	795 796 797	798 799 800	801 802 803	804 805 806	807 808 809	810 811 812	813 814 815	816 817 818	819 820 821	822 823 824	825 826 827	828 829 830	831 832 833	834 835 836	837 838 839	840 841 842	843 844 845	846 847 848	849 850 851	852 853 854	855 856 857	858 859 860	861 862 863	864 865 866	867 868 869	870 871 872	873 874 875	876 877 878	879 880 881	882 883 884	885 886 887	888 889 890	891 892 893	894 895 896	897 898 899	900 901 902	903 904 905	906 907 908	909 910 911	912 913 914	915 916 917	918 919 920	921 922 923	924 925 926	927 928 929	930 931 932	933 934 935	936 937 938	939 940 941	942 943 944	945 946 947	948 949 950	951 952 953	954 955 956	957 958 959	960 961 962	963 964 965	966 967 968	969 970 971	972 973 974	975 976 977	978 979 980	981 982 983	984 985 986	987 988 989	990 991 992	993 994 995	996 997 998	999 1000 1001	1002 1003 1004	1005 1006 1007	1008 1009 1010	1011 1012 1013	1014 1015 1016	1017 1018 1019	1020 1021 1022	1023 1024 1025	1026 1027 1028	1029 1030 1031	1032 1033 1034	1035 1036 1037	1038 1039 1040	1041 1042 1043	1044 1045 1046	1047 1048 1049	1050 1051 1052	1053 1054 1055	1056 1057 1058	1059 1060 1061	1062 1063 1064	1065 1066 1067	1068 1069 1070	1071 1072 1073	1074 1075 1076	1077 1078 1079	1080 1081 1082	1083 1084 1085	1086 1087 1088	1089 1090 1091	1092 1093 1094	1095 1096 1097	1098 1099 1100	1101 1102 1103	1104 1105 1106	1107 1108 1109	1110 1111 1112	1113 1114 1115	1116 1117 1118	1119 1120 1121	1122 1123 1124	1125 1126 1127	1128 1129 1130	1131 1132 1133	1134 1135 1136	1137 1138 1139	1140 1141 1142	1143 1144 1145	1146 1147 1148	1149 1150 1151	1152 1153 1154	1155 1156 1157	1158 1159 1160	1161 1162 1163	1164 1165 1166	1167 1168 1169	1170 1171 1172	1173 1174 1175	1176 1177 1178	1179 1180 1181	1182 1183 1184	1185 1186 1187	1188 1189 1190	1191 1192 1193	1194 1195 1196	1197 1198 1199	1200 1201 1202	1203 1204 1205	1206 1207 1208	1209 1210 1211	1212 1213 1214	1215 1216 1217	1218 1219 1220	1221 1222 1223	1224 1225 1226	1227 1228 1229	1230 1231 1232	1233 1234 1235	1236 1237 1238	1239 1240 1241	1242 1243 1244	1245 1246 1247	1248 1249 1250	1251 1252 1253	1254 1255 1256	1257 1258 1259	1260 1261 1262	1263 1264 1265	1266 1267 1268	1269 1270 1271	1272 1273 1274	1275 1276 1277	1278 1279 1280	1281 1282 1283	1284 1285 1286	1287 1288 1289	1290 1291 1292	1293 1294 1295	1296 1297 1298	1299 1300 1301	1302 1303 1304	1305 1306 1307	1308 1309 1310	1311 1312 1313	1314 1315 1316	1317 1318 1319	1320 1321 1322	1323 1324 1325	1326 1327 1328	1329 1330 1331	1332 1333 1334	1335 1336 1337	1338 1339 1340	1341 1342 1343	1344 1345 1346	1347 1348 1349	1350 1351 1352	1353 1354 1355	1356 1357 1358	1359 1360 1361	1362 1363 1364	1365 1366 1367	1368 1369 1370	1371 1372 1373	1374 1375 1376	1377 1378 1379	1380 1381 1382	1383 1384 1385	1386 1387 1388	1389 1390 1391	1392 1393 1394	1395 1396 1397	1398 1399 1400	1401 1402 1403	1404 1405 1406	1407 1408 1409	1410 1411 1412	1413 1414 1415	1416 1417 1418	1419 1420 1421	1422 1423 1424	1425 1426 1427	1428 1429 1430	1431 1432 1433	1434 1435 1436	1437 1438 1439	1440 1441 1442	1443 1444 1445	1446 1447 1448	1449 1450 1451	1452 1453 1454	1455 1456 1457	1458 1459 1460	1461 1462 1463	1464 1465 1466	1467 1468 1469	1470 1471 1472	1473 1474 1475	1476 1477 1478	1479 1480 1481	1482 1483 1484	1485 1486 1487	1488 1489 1490	1491 1492 1493	1494 1495 1496	1497 1498 1499	1500 1501 1502	1503 1504 1505	1506 1507 1508	1509 1510 1511	1512 1513 1514	1515 1516 1517	1518 1519 1520	1521 1522 1523	1524 1525 1526	1527 1528 1529	1530 1531 1532	1533 1534 1535	1536 1537 1538	1539 1540 1541	1542 1543 1544	1545 1546 1547	1548 1549 1550	1551 1552 1553	1554 1555 1556	1557 1558 1559	1560 1561 1562	1563 1564 1565	1566 1567 1568	1569 1570 1571	1572 1573 1574	1575 1576 1577	1578 1579 1580	1581 1582 1583	1584 1585 1586	1587 1588 1589	1590 1591 1592	1593 1594 1595	1596 1597 1598	1599 1600 1601	1602 1603 1604	1605 1606 1607	1608 1609 1610	1611 1612 1613	1614 1615 1616	1617 1618 1619	1620 1621 1622	1623 1624 1625	1626 1627 1628	1629 1630 1631	1632 1633 1634	1635 1636 1637	1638 1639 1640	1641 1642 1643	1644 1645 1646	1647 1648 1649	1650 1651 1652	1653 1654 1655	1656 1657 1658	1659 1660 1661	1662 1663 1664	1665 1666 1667	1668 1669 1670	1671 1672 1673	1674 1675 1676	1677 1678 1679	1680 1681 1682	1683 1684 1685	1686 1687 1688	1689 1690 1691	1692 1693 1694	1695 1696 1697	1698 1699 1700	1701 1702 1703	1704 1705 1706	1707 1708 1709	1710 1711 1712	1713 1714 1715	1716 1717 1718	1719 1720 1721	1722 1723 1724	1725 1726 1727	1728 1729 1730	1731 1732 1733	1734 1735 1736	1737 1738 1739	1740 1741 1742	1743 1744 1745	1746 1747 1748	1749 1750 1751	1752 1753 1754	1755 1756 1757	1758 1759 1760	1761 1762 1763	1764 1765 1766	1767 1768 1769	1770 1771 1772	1773 1774 1775	1776 1777 1778	1779 1780 1781	1782 1783 1784	1785 1786 1787	1788 1789 1790	1791 1792 1793	1794 1795 1796	1797 1798 1799	1800 1801 1802	1803 1804 1805	1806 1807 1808	1809 1810 1811	1812 1813 1814	1815 1816 1817	1818 1819 1820	1821 1822 1823	1824 1825 1826	1827 1828 1829	1830 1831 1832	1833 1834 1835	1836 1837 1838	1839 1840 1841	1842 1843 1844	1845 1846 1847	1848 1849 1850	1851 1852 1853	1854 1855 1856	1857 1858 1859	1860 1861 1862	1863 1864 1865	1866 1867 1868	1869 1870 1871	1872 1873 1874	1875 1876 1877	1878 1879 1880	1881 1882 1883	1884 1885 1886	1887 1888 1889	1890 1891 1892	1893 1894 1895	1896 1897 1898	1899 1900 1901	1902 1903 1904	1905 1906 1907	1908 1909 1910	1911 1912 1913	1914 1915 1916	1917 1918 1919	1920 1921 1922	1923 1924 1925	1926 1927 1928	1929 1930 1931	1932 1933 1934	1935 1936 1937	1938 1939 1940	1941 1942 1943	1944 1945 1946	1947 1948 1949	1950 1951 1952	1953 1954 1955	1956 1957 1958	1959 1960 1961	1962 1963 1964	1965 1966 1967	1968 1969 1970	1971 1972 1973	1974 1975 1976	1977 1978 1979	1980 1981 1982	1983 1984 1985	1986 1987 1988	1989 1990 1991	1992 1993 1994	1995 1996 1997	1998 1999 2000	2001 2002 2003	2004 2005 2006	2007 2008 2009	2010 2011 2012	2013 2014 2015	2016 2017 2018	2019 2020 2021	2022 2023 2024	2025 2026 2027	2028 2029 2030	2031 2032 2033	2034 2035 2036	2037 2038 2039	2040 2041 2042	2043 2044 2045	2046 2047 2048	2049 2050 2051	2052 2053 2054	2055 2056 2057	2058 2059 2060	2061 2062 2063	2064 2065 2066	2067 2068 2069	2070 2071 2072	2073 2074 2075	2076 2077 2078	2079 2080 2081	2082 2083 2084	2085 2086 2087	2088 2089 2090	2091 2092 2093	2094 2095 2096	2097 2098 2099	2100 2101 2102	2103 2104 2105	2106 2107 2108	2109 2110 2111	2112 2113 2114	2115 2116 2117	2118 2119 2120	2121 2122 2123	2124 2125 2126	2127 2128 2129	2130 2131 2132	2133 2134 2135	2136 2137 2138	2139 2140 2141	2142 2143 2144	2145 2146 2147	2148 2149 2150	2151 2152 2153	2154 2155 2156	2157 2158 2159	2160 2161 2162	2163 2164 2165	2166 2167 2168	2169 2170 2171	2172 2173 2174	2175 2176 2177	2178 2179 2180	2181 2182 2183	2184 2185 2186	2187 2188 2189	2190 2191 2192	2193 2194 2195	2196 2197 2198	2199 2200 2201	2202 2203 2204	2205 2206 2207	2208 2209 2210	2211 2212 2213	2214 2215 2216	2217 2218 2219	2220 2221 2222	2223 2224 2225	2226 2227 2228	2229 2230 2231	2232 2233 2234	2235 2236 2237	2238 2239 2240	2241 2242 2243	2244 2245 2246	2247 2248 2249	2250 2251 2252	2253 2254 2255	2256 2257 2258	2259 2260 2261	2262 2263 2264	2265 2266 2267	2268 2269 2270	2271 2272 2273	2274 2275 2276	2277 2278 2279	2280 2281 2282	2283 2284 2285	2286 2287 2288	2289 2290 2291	2292 2293 2294	2295 2296 2297	2298 2299 2300	2301 2302 2303	2304 2305 2306	2307 2308 2309	2310 2311 2312	2313 2314 2315	2316 2317 2318	2319 2320 2321	2322 2323 2324	2325 2326 2327	2328 2329 2330	2331 2332 2333	2334 2335 2336	2337 2338 2339	2340 2341 2342	2343 2344 2345	2346 2347 2348	2349 2350 2351	2352 2353 2354	2355 2356 2357	2358 2359 2360	2361 2362 2363	2364 2365 2366	2367 2368 2369	2370 2371 2372	2373 2374 2375	2376 2377 2378	2379 2380 2381	2382 2383 2384	2385 2386 2387	2388 2389 2390	2391 2392 2393	2394 2395 2396	2397 2398 2399	2400 2401 2402	2403 2404 2405	2406 2407 2408	2409 2410 2411	2412 2413 2414	2415 2416 2417	2418 2419 2420	2421 2422 2423	2424 2425 2426	2427 2428 2429	2430 2431 2432	2433 2434 2435	2436 2437 2438	2439 2440 2441	2442 2443 2444	2445 2446 2447	2448 2449 2450	2451 2452 2453	2454 2455 2456	2457 2458 2459	2460 2461 2462	2463 2464 2465	2466 2467 2468	2469 2470 2471	2472 2473 2474	2475 2476 2477	2478 2479 2480	2481 2482 2483	2484 2485 2486	2487 2488 2489	2490 2491 2492	2493 2494 2495	2496 2497 2498	2499 2500 2501	2502 2503 2504	2505 2506 2507	2508 2509 2510	2511 2512 2513	2514 2515 2516	2517 2518 2519	2520 2521 2522	2523 2524 2525	2526 2527 2528	2529 2530 2531	2532 2533 2534	2535 2536 2537	2538 2539 2540	2541 2542 2543	2544 2545 2546	2547 2548 2549	2550 2551 2552	2553 2554 2555	2556 2557 2558	2559 2560 2561	2562 2563 2564	2565 2566 2567	2568 2569 2570	2571 2572 2573	2574 2575 2576	2577 2578 2579	2580 2581 2582	2583 2584 2585	2586 2587 2588	2589 2590 2591	2592 2593 2594	2595 2596 2597	2598 2599 2600	2601 2602 2603	2604 2605 2606	2607 2608 2609	2610 2611 2612	2613 2614 2615	2616 2617 2618	2619 2620 2621	2622 2623 2624	2625 2626 2627	2628 2629 2630	2631 2632 2633	2634 2635 2636	2637 2638 2639	2640 2641 2642	2643 2644 2645	2646 2647 2648	2649 2650 2651	2652 2653 2654	2655 2656 2657	2658 2659 2660	2661 2662 2663	2664 2665 2666	2667 2668 2669	2670 2671 2672	2673 2674 2675	2676 2677 2678	2679 2680 2681	2682 2683 2684	2685 2686 2687	2688 2689 2690	2691 2692 2693	2694 2695 2696	2697 2698 2699	2700 2701 2702	2703 2704 2705	2706 2707 2708	2709 2710 2711	2712 2713 2714	2715 2716 2717	2718 2719 2720	2721 2722 2723	2724 2725 2726	2727 2728 2729	2730 2731 2732	2733 2734 2735	2736 2737 2738	2739 2740 2741	2742 2743 2744	2745 2746 2747	2748 2749 2750	2751 2752 2753	2754 2755 2756	2757 2758 2759	2760 2761 2762	2763 2764 2765	2766 2767 2768	2769 2770 2771	2772 2773 2774	2775 2776 2777	2778 2779 2780	2781 2782 2783	2784 2785 2786	2787 2788 2789	2790 2791 2792	2793 2794 2795	2796 2797 2798	2799 2800 2801	2802 2803 2804	2805 2806 2807	2808 2809 2810	2811 2812 2813	2814 2815 2816	2817 2818 2819	2820 2821 2822	2823 2824 2825	2826 2827 2828	2829 2830 2831	2832 2833 2834	2835 2836 2837	2838 2839 2840	2841 2842 2843	2844 2845 2846	2847 2848 2849	2850 2851 2852	2853 2854 2855	2856 2857 2858	2859 2860 2861	2862 2863 2864	2865 2866 2867	2868 2869 2870	2871 2872 2873	2874 2875 2876	2877 2878 2879]\n" to:out_file
		format "\"point P\" [0.0649438 -0.121756 0.160653	0.0649438 -0.105966 0.159098	0.0618633 -0.10627 0.159098	0.0649438 -0.121756 0.160653	0.0618633 -0.10627 0.159098	0.0589013 -0.107168 0.159098	0.0649438 -0.121756 0.160653	0.0589013 -0.107168 0.159098	0.0561714 -0.108627 0.159098	0.0649438 -0.121756 0.160653	0.0561714 -0.108627 0.159098	0.0537786 -0.110591 0.159098	0.0649438 -0.121756 0.160653	0.0537786 -0.110591 0.159098	0.0518149 -0.112984 0.159098	0.0649438 -0.121756 0.160653	0.0518149 -0.112984 0.159098	0.0503558 -0.115714 0.159098	0.0649438 -0.121756 0.160653	0.0503558 -0.115714 0.159098	0.0494573 -0.118676 0.159098	0.0649438 -0.121756 0.160653	0.0494573 -0.118676 0.159098	0.0491539 -0.121756 0.159098	0.0649438 -0.121756 0.160653	0.0491539 -0.121756 0.159098	0.0494573 -0.124837 0.159098	0.0649438 -0.121756 0.160653	0.0494573 -0.124837 0.159098	0.0503558 -0.127799 0.159098	0.0649438 -0.121756 0.160653	0.0503558 -0.127799 0.159098	0.051815 -0.130528 0.159098	0.0649438 -0.121756 0.160653	0.051815 -0.130528 0.159098	0.0537786 -0.132921 0.159098	0.0649438 -0.121756 0.160653	0.0537786 -0.132921 0.159098	0.0561714 -0.134885 0.159098	0.0649438 -0.121756 0.160653	0.0561714 -0.134885 0.159098	0.0589013 -0.136344 0.159098	0.0649438 -0.121756 0.160653	0.0589013 -0.136344 0.159098	0.0618634 -0.137243 0.159098	0.0649438 -0.121756 0.160653	0.0618634 -0.137243 0.159098	0.0649438 -0.137546 0.159098	0.0649438 -0.121756 0.160653	0.0649438 -0.137546 0.159098	0.0680243 -0.137243 0.159098	0.0649438 -0.121756 0.160653	0.0680243 -0.137243 0.159098	0.0709864 -0.136344 0.159098	0.0649438 -0.121756 0.160653	0.0709864 -0.136344 0.159098	0.0737163 -0.134885 0.159098	0.0649438 -0.121756 0.160653	0.0737163 -0.134885 0.159098	0.076109 -0.132921 0.159098	0.0649438 -0.121756 0.160653	0.076109 -0.132921 0.159098	0.0780727 -0.130528 0.159098	0.0649438 -0.121756 0.160653	0.0780727 -0.130528 0.159098	0.0795318 -0.127799 0.159098	0.0649438 -0.121756 0.160653	0.0795318 -0.127799 0.159098	0.0804304 -0.124837 0.159098	0.0649438 -0.121756 0.160653	0.0804304 -0.124837 0.159098	0.0807338 -0.121756 0.159098	0.0649438 -0.121756 0.160653	0.0807338 -0.121756 0.159098	0.0804304 -0.118676 0.159098	0.0649438 -0.121756 0.160653	0.0804304 -0.118676 0.159098	0.0795318 -0.115713 0.159098	0.0649438 -0.121756 0.160653	0.0795318 -0.115713 0.159098	0.0780727 -0.112984 0.159098	0.0649438 -0.121756 0.160653	0.0780727 -0.112984 0.159098	0.076109 -0.110591 0.159098	0.0649438 -0.121756 0.160653	0.076109 -0.110591 0.159098	0.0737162 -0.108627 0.159098	0.0649438 -0.121756 0.160653	0.0737162 -0.108627 0.159098	0.0709863 -0.107168 0.159098	0.0649438 -0.121756 0.160653	0.0709863 -0.107168 0.159098	0.0680242 -0.10627 0.159098	0.0649438 -0.121756 0.160653	0.0680242 -0.10627 0.159098	0.0649438 -0.105966 0.159098	0.0649438 -0.105966 0.159098	0.0649438 -0.090783 0.154492	0.0589013 -0.0913781 0.154492	0.0649438 -0.105966 0.159098	0.0589013 -0.0913781 0.154492	0.0618633 -0.10627 0.159098	0.0618633 -0.10627 0.159098	0.0589013 -0.0913781 0.154492	0.0530909 -0.0931406 0.154492	0.0618633 -0.10627 0.159098	0.0530909 -0.0931406 0.154492	0.0589013 -0.107168 0.159098	0.0589013 -0.107168 0.159098	0.0530909 -0.0931406 0.154492	0.0477361 -0.0960029 0.154492	0.0589013 -0.107168 0.159098	0.0477361 -0.0960029 0.154492	0.0561714 -0.108627 0.159098	0.0561714 -0.108627 0.159098	0.0477361 -0.0960029 0.154492	0.0430425 -0.0998548 0.154492	0.0561714 -0.108627 0.159098	0.0430425 -0.0998548 0.154492	0.0537786 -0.110591 0.159098	0.0537786 -0.110591 0.159098	0.0430425 -0.0998548 0.154492	0.0391906 -0.104548 0.154492	0.0537786 -0.110591 0.159098	0.0391906 -0.104548 0.154492	0.0518149 -0.112984 0.159098	0.0518149 -0.112984 0.159098	0.0391906 -0.104548 0.154492	0.0363284 -0.109903 0.154492	0.0518149 -0.112984 0.159098	0.0363284 -0.109903 0.154492	0.0503558 -0.115714 0.159098	0.0503558 -0.115714 0.159098	0.0363284 -0.109903 0.154492	0.0345659 -0.115714 0.154492	0.0503558 -0.115714 0.159098	0.0345659 -0.115714 0.154492	0.0494573 -0.118676 0.159098	0.0494573 -0.118676 0.159098	0.0345659 -0.115714 0.154492	0.0339707 -0.121756 0.154492	0.0494573 -0.118676 0.159098	0.0339707 -0.121756 0.154492	0.0491539 -0.121756 0.159098	0.0491539 -0.121756 0.159098	0.0339707 -0.121756 0.154492	0.0345659 -0.127799 0.154492	0.0491539 -0.121756 0.159098	0.0345659 -0.127799 0.154492	0.0494573 -0.124837 0.159098	0.0494573 -0.124837 0.159098	0.0345659 -0.127799 0.154492	0.0363284 -0.133609 0.154492	0.0494573 -0.124837 0.159098	0.0363284 -0.133609 0.154492	0.0503558 -0.127799 0.159098	0.0503558 -0.127799 0.159098	0.0363284 -0.133609 0.154492	0.0391906 -0.138964 0.154492	0.0503558 -0.127799 0.159098	0.0391906 -0.138964 0.154492	0.051815 -0.130528 0.159098	0.051815 -0.130528 0.159098	0.0391906 -0.138964 0.154492	0.0430425 -0.143657 0.154492	0.051815 -0.130528 0.159098	0.0430425 -0.143657 0.154492	0.0537786 -0.132921 0.159098	0.0537786 -0.132921 0.159098	0.0430425 -0.143657 0.154492	0.0477361 -0.147509 0.154492	0.0537786 -0.132921 0.159098	0.0477361 -0.147509 0.154492	0.0561714 -0.134885 0.159098	0.0561714 -0.134885 0.159098	0.0477361 -0.147509 0.154492	0.0530909 -0.150371 0.154492	0.0561714 -0.134885 0.159098	0.0530909 -0.150371 0.154492	0.0589013 -0.136344 0.159098	0.0589013 -0.136344 0.159098	0.0530909 -0.150371 0.154492	0.0589013 -0.152134 0.154492	0.0589013 -0.136344 0.159098	0.0589013 -0.152134 0.154492	0.0618634 -0.137243 0.159098	0.0618634 -0.137243 0.159098	0.0589013 -0.152134 0.154492	0.0649438 -0.152729 0.154492	0.0618634 -0.137243 0.159098	0.0649438 -0.152729 0.154492	0.0649438 -0.137546 0.159098	0.0649438 -0.137546 0.159098	0.0649438 -0.152729 0.154492	0.0709864 -0.152134 0.154492	0.0649438 -0.137546 0.159098	0.0709864 -0.152134 0.154492	0.0680243 -0.137243 0.159098	0.0680243 -0.137243 0.159098	0.0709864 -0.152134 0.154492	0.0767967 -0.150371 0.154492	0.0680243 -0.137243 0.159098	0.0767967 -0.150371 0.154492	0.0709864 -0.136344 0.159098	0.0709864 -0.136344 0.159098	0.0767967 -0.150371 0.154492	0.0821516 -0.147509 0.154492	0.0709864 -0.136344 0.159098	0.0821516 -0.147509 0.154492	0.0737163 -0.134885 0.159098	0.0737163 -0.134885 0.159098	0.0821516 -0.147509 0.154492	0.0868451 -0.143657 0.154492	0.0737163 -0.134885 0.159098	0.0868451 -0.143657 0.154492	0.076109 -0.132921 0.159098	0.076109 -0.132921 0.159098	0.0868451 -0.143657 0.154492	0.090697 -0.138964 0.154492	0.076109 -0.132921 0.159098	0.090697 -0.138964 0.154492	0.0780727 -0.130528 0.159098	0.0780727 -0.130528 0.159098	0.090697 -0.138964 0.154492	0.0935592 -0.133609 0.154492	0.0780727 -0.130528 0.159098	0.0935592 -0.133609 0.154492	0.0795318 -0.127799 0.159098	0.0795318 -0.127799 0.159098	0.0935592 -0.133609 0.154492	0.0953218 -0.127799 0.154492	0.0795318 -0.127799 0.159098	0.0953218 -0.127799 0.154492	0.0804304 -0.124837 0.159098	0.0804304 -0.124837 0.159098	0.0953218 -0.127799 0.154492	0.0959169 -0.121756 0.154492	0.0804304 -0.124837 0.159098	0.0959169 -0.121756 0.154492	0.0807338 -0.121756 0.159098	0.0807338 -0.121756 0.159098	0.0959169 -0.121756 0.154492	0.0953218 -0.115713 0.154492	0.0807338 -0.121756 0.159098	0.0953218 -0.115713 0.154492	0.0804304 -0.118676 0.159098	0.0804304 -0.118676 0.159098	0.0953218 -0.115713 0.154492	0.0935592 -0.109903 0.154492	0.0804304 -0.118676 0.159098	0.0935592 -0.109903 0.154492	0.0795318 -0.115713 0.159098	0.0795318 -0.115713 0.159098	0.0935592 -0.109903 0.154492	0.090697 -0.104548 0.154492	0.0795318 -0.115713 0.159098	0.090697 -0.104548 0.154492	0.0780727 -0.112984 0.159098	0.0780727 -0.112984 0.159098	0.090697 -0.104548 0.154492	0.0868451 -0.0998547 0.154492	0.0780727 -0.112984 0.159098	0.0868451 -0.0998547 0.154492	0.076109 -0.110591 0.159098	0.076109 -0.110591 0.159098	0.0868451 -0.0998547 0.154492	0.0821515 -0.0960028 0.154492	0.076109 -0.110591 0.159098	0.0821515 -0.0960028 0.154492	0.0737162 -0.108627 0.159098	0.0737162 -0.108627 0.159098	0.0821515 -0.0960028 0.154492	0.0767966 -0.0931406 0.154492	0.0737162 -0.108627 0.159098	0.0767966 -0.0931406 0.154492	0.0709863 -0.107168 0.159098	0.0709863 -0.107168 0.159098	0.0767966 -0.0931406 0.154492	0.0709863 -0.0913781 0.154492	0.0709863 -0.107168 0.159098	0.0709863 -0.0913781 0.154492	0.0680242 -0.10627 0.159098	0.0680242 -0.10627 0.159098	0.0709863 -0.0913781 0.154492	0.0649438 -0.090783 0.154492	0.0680242 -0.10627 0.159098	0.0649438 -0.090783 0.154492	0.0649438 -0.105966 0.159098	0.0649438 -0.090783 0.154492	0.0649438 -0.0767901 0.147012	0.0561714 -0.0776541 0.147012	0.0649438 -0.090783 0.154492	0.0561714 -0.0776541 0.147012	0.0589013 -0.0913781 0.154492	0.0589013 -0.0913781 0.154492	0.0561714 -0.0776541 0.147012	0.0477361 -0.0802129 0.147012	0.0589013 -0.0913781 0.154492	0.0477361 -0.0802129 0.147012	0.0530909 -0.0931406 0.154492	0.0530909 -0.0931406 0.154492	0.0477361 -0.0802129 0.147012	0.0399621 -0.0843682 0.147012	0.0530909 -0.0931406 0.154492	0.0399621 -0.0843682 0.147012	0.0477361 -0.0960029 0.154492	0.0477361 -0.0960029 0.154492	0.0399621 -0.0843682 0.147012	0.0331481 -0.0899603 0.147012	0.0477361 -0.0960029 0.154492	0.0331481 -0.0899603 0.147012	0.0430425 -0.0998548 0.154492	0.0430425 -0.0998548 0.154492	0.0331481 -0.0899603 0.147012	0.027556 -0.0967743 0.147012	0.0430425 -0.0998548 0.154492	0.027556 -0.0967743 0.147012	0.0391906 -0.104548 0.154492	0.0391906 -0.104548 0.154492	0.027556 -0.0967743 0.147012	0.0234007 -0.104548 0.147012	0.0391906 -0.104548 0.154492	0.0234007 -0.104548 0.147012	0.0363284 -0.109903 0.154492	0.0363284 -0.109903 0.154492	0.0234007 -0.104548 0.147012	0.0208419 -0.112984 0.147012	0.0363284 -0.109903 0.154492	0.0208419 -0.112984 0.147012	0.0345659 -0.115714 0.154492	0.0345659 -0.115714 0.154492	0.0208419 -0.112984 0.147012	0.0199778 -0.121756 0.147012	0.0345659 -0.115714 0.154492	0.0199778 -0.121756 0.147012	0.0339707 -0.121756 0.154492	0.0339707 -0.121756 0.154492	0.0199778 -0.121756 0.147012	0.0208419 -0.130529 0.147012	0.0339707 -0.121756 0.154492	0.0208419 -0.130529 0.147012	0.0345659 -0.127799 0.154492	0.0345659 -0.127799 0.154492	0.0208419 -0.130529 0.147012	0.0234007 -0.138964 0.147012	0.0345659 -0.127799 0.154492	0.0234007 -0.138964 0.147012	0.0363284 -0.133609 0.154492	0.0363284 -0.133609 0.154492	0.0234007 -0.138964 0.147012	0.027556 -0.146738 0.147012	0.0363284 -0.133609 0.154492	0.027556 -0.146738 0.147012	0.0391906 -0.138964 0.154492	0.0391906 -0.138964 0.154492	0.027556 -0.146738 0.147012	0.0331481 -0.153552 0.147012	0.0391906 -0.138964 0.154492	0.0331481 -0.153552 0.147012	0.0430425 -0.143657 0.154492	0.0430425 -0.143657 0.154492	0.0331481 -0.153552 0.147012	0.0399621 -0.159144 0.147012	0.0430425 -0.143657 0.154492	0.0399621 -0.159144 0.147012	0.0477361 -0.147509 0.154492	0.0477361 -0.147509 0.154492	0.0399621 -0.159144 0.147012	0.0477361 -0.163299 0.147012	0.0477361 -0.147509 0.154492	0.0477361 -0.163299 0.147012	0.0530909 -0.150371 0.154492	0.0530909 -0.150371 0.154492	0.0477361 -0.163299 0.147012	0.0561714 -0.165858 0.147012	0.0530909 -0.150371 0.154492	0.0561714 -0.165858 0.147012	0.0589013 -0.152134 0.154492	0.0589013 -0.152134 0.154492	0.0561714 -0.165858 0.147012	0.0649439 -0.166722 0.147012	0.0589013 -0.152134 0.154492	0.0649439 -0.166722 0.147012	0.0649438 -0.152729 0.154492	0.0649438 -0.152729 0.154492	0.0649439 -0.166722 0.147012	0.0737163 -0.165858 0.147012	0.0649438 -0.152729 0.154492	0.0737163 -0.165858 0.147012	0.0709864 -0.152134 0.154492	0.0709864 -0.152134 0.154492	0.0737163 -0.165858 0.147012	0.0821516 -0.163299 0.147012	0.0709864 -0.152134 0.154492	0.0821516 -0.163299 0.147012	0.0767967 -0.150371 0.154492	0.0767967 -0.150371 0.154492	0.0821516 -0.163299 0.147012	0.0899256 -0.159144 0.147012	0.0767967 -0.150371 0.154492	0.0899256 -0.159144 0.147012	0.0821516 -0.147509 0.154492	0.0821516 -0.147509 0.154492	0.0899256 -0.159144 0.147012	0.0967396 -0.153552 0.147012	0.0821516 -0.147509 0.154492	0.0967396 -0.153552 0.147012	0.0868451 -0.143657 0.154492	0.0868451 -0.143657 0.154492	0.0967396 -0.153552 0.147012	0.102332 -0.146738 0.147012	0.0868451 -0.143657 0.154492	0.102332 -0.146738 0.147012	0.090697 -0.138964 0.154492	0.090697 -0.138964 0.154492	0.102332 -0.146738 0.147012	0.106487 -0.138964 0.147012	0.090697 -0.138964 0.154492	0.106487 -0.138964 0.147012	0.0935592 -0.133609 0.154492	0.0935592 -0.133609 0.154492	0.106487 -0.138964 0.147012	0.109046 -0.130528 0.147012	0.0935592 -0.133609 0.154492	0.109046 -0.130528 0.147012	0.0953218 -0.127799 0.154492	0.0953218 -0.127799 0.154492	0.109046 -0.130528 0.147012	0.10991 -0.121756 0.147012	0.0953218 -0.127799 0.154492	0.10991 -0.121756 0.147012	0.0959169 -0.121756 0.154492	0.0959169 -0.121756 0.154492	0.10991 -0.121756 0.147012	0.109046 -0.112984 0.147012	0.0959169 -0.121756 0.154492	0.109046 -0.112984 0.147012	0.0953218 -0.115713 0.154492	0.0953218 -0.115713 0.154492	0.109046 -0.112984 0.147012	0.106487 -0.104548 0.147012	0.0953218 -0.115713 0.154492	0.106487 -0.104548 0.147012	0.0935592 -0.109903 0.154492	0.0935592 -0.109903 0.154492	0.106487 -0.104548 0.147012	0.102332 -0.0967742 0.147012	0.0935592 -0.109903 0.154492	0.102332 -0.0967742 0.147012	0.090697 -0.104548 0.154492	0.090697 -0.104548 0.154492	0.102332 -0.0967742 0.147012	0.0967395 -0.0899603 0.147012	0.090697 -0.104548 0.154492	0.0967395 -0.0899603 0.147012	0.0868451 -0.0998547 0.154492	0.0868451 -0.0998547 0.154492	0.0967395 -0.0899603 0.147012	0.0899255 -0.0843682 0.147012	0.0868451 -0.0998547 0.154492	0.0899255 -0.0843682 0.147012	0.0821515 -0.0960028 0.154492	0.0821515 -0.0960028 0.154492	0.0899255 -0.0843682 0.147012	0.0821515 -0.0802129 0.147012	0.0821515 -0.0960028 0.154492	0.0821515 -0.0802129 0.147012	0.0767966 -0.0931406 0.154492	0.0767966 -0.0931406 0.154492	0.0821515 -0.0802129 0.147012	0.0737161 -0.0776541 0.147012	0.0767966 -0.0931406 0.154492	0.0737161 -0.0776541 0.147012	0.0709863 -0.0913781 0.154492	0.0709863 -0.0913781 0.154492	0.0737161 -0.0776541 0.147012	0.0649438 -0.0767901 0.147012	0.0709863 -0.0913781 0.154492	0.0649438 -0.0767901 0.147012	0.0649438 -0.090783 0.154492	0.0649438 -0.0767901 0.147012	0.0649438 -0.0645252 0.136947	0.0537786 -0.0656249 0.136947	0.0649438 -0.0767901 0.147012	0.0537786 -0.0656249 0.136947	0.0561714 -0.0776541 0.147012	0.0561714 -0.0776541 0.147012	0.0537786 -0.0656249 0.136947	0.0430425 -0.0688817 0.136947	0.0561714 -0.0776541 0.147012	0.0430425 -0.0688817 0.136947	0.0477361 -0.0802129 0.147012	0.0477361 -0.0802129 0.147012	0.0430425 -0.0688817 0.136947	0.0331481 -0.0741704 0.136947	0.0477361 -0.0802129 0.147012	0.0331481 -0.0741704 0.136947	0.0399621 -0.0843682 0.147012	0.0399621 -0.0843682 0.147012	0.0331481 -0.0741704 0.136947	0.0244755 -0.0812878 0.136947	0.0399621 -0.0843682 0.147012	0.0244755 -0.0812878 0.136947	0.0331481 -0.0899603 0.147012	0.0331481 -0.0899603 0.147012	0.0244755 -0.0812878 0.136947	0.0173581 -0.0899603 0.136947	0.0331481 -0.0899603 0.147012	0.0173581 -0.0899603 0.136947	0.027556 -0.0967743 0.147012	0.027556 -0.0967743 0.147012	0.0173581 -0.0899603 0.136947	0.0120694 -0.0998548 0.136947	0.027556 -0.0967743 0.147012	0.0120694 -0.0998548 0.136947	0.0234007 -0.104548 0.147012	0.0234007 -0.104548 0.147012	0.0120694 -0.0998548 0.136947	0.00881266 -0.110591 0.136947	0.0234007 -0.104548 0.147012	0.00881266 -0.110591 0.136947	0.0208419 -0.112984 0.147012	0.0208419 -0.112984 0.147012	0.00881266 -0.110591 0.136947	0.00771299 -0.121756 0.136947	0.0208419 -0.112984 0.147012	0.00771299 -0.121756 0.136947	0.0199778 -0.121756 0.147012	0.0199778 -0.121756 0.147012	0.00771299 -0.121756 0.136947	0.00881267 -0.132921 0.136947	0.0199778 -0.121756 0.147012	0.00881267 -0.132921 0.136947	0.0208419 -0.130529 0.147012	0.0208419 -0.130529 0.147012	0.00881267 -0.132921 0.136947	0.0120694 -0.143657 0.136947	0.0208419 -0.130529 0.147012	0.0120694 -0.143657 0.136947	0.0234007 -0.138964 0.147012	0.0234007 -0.138964 0.147012	0.0120694 -0.143657 0.136947	0.0173581 -0.153552 0.136947	0.0234007 -0.138964 0.147012	0.0173581 -0.153552 0.136947	0.027556 -0.146738 0.147012	0.027556 -0.146738 0.147012	0.0173581 -0.153552 0.136947	0.0244755 -0.162224 0.136947	0.027556 -0.146738 0.147012	0.0244755 -0.162224 0.136947	0.0331481 -0.153552 0.147012	0.0331481 -0.153552 0.147012	0.0244755 -0.162224 0.136947	0.0331481 -0.169342 0.136947	0.0331481 -0.153552 0.147012	0.0331481 -0.169342 0.136947	0.0399621 -0.159144 0.147012	0.0399621 -0.159144 0.147012	0.0331481 -0.169342 0.136947	0.0430426 -0.17463 0.136947	0.0399621 -0.159144 0.147012	0.0430426 -0.17463 0.136947	0.0477361 -0.163299 0.147012	0.0477361 -0.163299 0.147012	0.0430426 -0.17463 0.136947	0.0537787 -0.177887 0.136947	0.0477361 -0.163299 0.147012	0.0537787 -0.177887 0.136947	0.0561714 -0.165858 0.147012	0.0561714 -0.165858 0.147012	0.0537787 -0.177887 0.136947	0.0649439 -0.178987 0.136947	0.0561714 -0.165858 0.147012	0.0649439 -0.178987 0.136947	0.0649439 -0.166722 0.147012	0.0649439 -0.166722 0.147012	0.0649439 -0.178987 0.136947	0.0761091 -0.177887 0.136947	0.0649439 -0.166722 0.147012	0.0761091 -0.177887 0.136947	0.0737163 -0.165858 0.147012	0.0737163 -0.165858 0.147012	0.0761091 -0.177887 0.136947	0.0868452 -0.17463 0.136947	0.0737163 -0.165858 0.147012	0.0868452 -0.17463 0.136947	0.0821516 -0.163299 0.147012	0.0821516 -0.163299 0.147012	0.0868452 -0.17463 0.136947	0.0967396 -0.169342 0.136947	0.0821516 -0.163299 0.147012	0.0967396 -0.169342 0.136947	0.0899256 -0.159144 0.147012	0.0899256 -0.159144 0.147012	0.0967396 -0.169342 0.136947	0.105412 -0.162224 0.136947	0.0899256 -0.159144 0.147012	0.105412 -0.162224 0.136947	0.0967396 -0.153552 0.147012	0.0967396 -0.153552 0.147012	0.105412 -0.162224 0.136947	0.11253 -0.153552 0.136947	0.0967396 -0.153552 0.147012	0.11253 -0.153552 0.136947	0.102332 -0.146738 0.147012	0.102332 -0.146738 0.147012	0.11253 -0.153552 0.136947	0.117818 -0.143657 0.136947	0.102332 -0.146738 0.147012	0.117818 -0.143657 0.136947	0.106487 -0.138964 0.147012	0.106487 -0.138964 0.147012	0.117818 -0.143657 0.136947	0.121075 -0.132921 0.136947	0.106487 -0.138964 0.147012	0.121075 -0.132921 0.136947	0.109046 -0.130528 0.147012	0.109046 -0.130528 0.147012	0.121075 -0.132921 0.136947	0.122175 -0.121756 0.136947	0.109046 -0.130528 0.147012	0.122175 -0.121756 0.136947	0.10991 -0.121756 0.147012	0.10991 -0.121756 0.147012	0.122175 -0.121756 0.136947	0.121075 -0.110591 0.136947	0.10991 -0.121756 0.147012	0.121075 -0.110591 0.136947	0.109046 -0.112984 0.147012	0.109046 -0.112984 0.147012	0.121075 -0.110591 0.136947	0.117818 -0.0998547 0.136947	0.109046 -0.112984 0.147012	0.117818 -0.0998547 0.136947	0.106487 -0.104548 0.147012	0.106487 -0.104548 0.147012	0.117818 -0.0998547 0.136947	0.112529 -0.0899602 0.136947	0.106487 -0.104548 0.147012	0.112529 -0.0899602 0.136947	0.102332 -0.0967742 0.147012	0.102332 -0.0967742 0.147012	0.112529 -0.0899602 0.136947	0.105412 -0.0812877 0.136947	0.102332 -0.0967742 0.147012	0.105412 -0.0812877 0.136947	0.0967395 -0.0899603 0.147012	0.0967395 -0.0899603 0.147012	0.105412 -0.0812877 0.136947	0.0967395 -0.0741703 0.136947	0.0967395 -0.0899603 0.147012	0.0967395 -0.0741703 0.136947	0.0899255 -0.0843682 0.147012	0.0899255 -0.0843682 0.147012	0.0967395 -0.0741703 0.136947	0.086845 -0.0688816 0.136947	0.0899255 -0.0843682 0.147012	0.086845 -0.0688816 0.136947	0.0821515 -0.0802129 0.147012	0.0821515 -0.0802129 0.147012	0.086845 -0.0688816 0.136947	0.0761089 -0.0656249 0.136947	0.0821515 -0.0802129 0.147012	0.0761089 -0.0656249 0.136947	0.0737161 -0.0776541 0.147012	0.0737161 -0.0776541 0.147012	0.0761089 -0.0656249 0.136947	0.0649438 -0.0645252 0.136947	0.0737161 -0.0776541 0.147012	0.0649438 -0.0645252 0.136947	0.0649438 -0.0767901 0.147012	0.0649438 -0.0645252 0.136947	0.0649438 -0.0544597 0.124682	0.051815 -0.0557528 0.124682	0.0649438 -0.0645252 0.136947	0.051815 -0.0557528 0.124682	0.0537786 -0.0656249 0.136947	0.0537786 -0.0656249 0.136947	0.051815 -0.0557528 0.124682	0.0391906 -0.0595824 0.124682	0.0537786 -0.0656249 0.136947	0.0391906 -0.0595824 0.124682	0.0430425 -0.0688817 0.136947	0.0430425 -0.0688817 0.136947	0.0391906 -0.0595824 0.124682	0.027556 -0.0658012 0.124682	0.0430425 -0.0688817 0.136947	0.027556 -0.0658012 0.124682	0.0331481 -0.0741704 0.136947	0.0331481 -0.0741704 0.136947	0.027556 -0.0658012 0.124682	0.0173581 -0.0741704 0.124682	0.0331481 -0.0741704 0.136947	0.0173581 -0.0741704 0.124682	0.0244755 -0.0812878 0.136947	0.0244755 -0.0812878 0.136947	0.0173581 -0.0741704 0.124682	0.00898895 -0.0843682 0.124682	0.0244755 -0.0812878 0.136947	0.00898895 -0.0843682 0.124682	0.0173581 -0.0899603 0.136947	0.0173581 -0.0899603 0.136947	0.00898895 -0.0843682 0.124682	0.0027701 -0.0960029 0.124682	0.0173581 -0.0899603 0.136947	0.0027701 -0.0960029 0.124682	0.0120694 -0.0998548 0.136947	0.0120694 -0.0998548 0.136947	0.0027701 -0.0960029 0.124682	-0.00105944 -0.108627 0.124682	0.0120694 -0.0998548 0.136947	-0.00105944 -0.108627 0.124682	0.00881266 -0.110591 0.136947	0.00881266 -0.110591 0.136947	-0.00105944 -0.108627 0.124682	-0.00235252 -0.121756 0.124682	0.00881266 -0.110591 0.136947	-0.00235252 -0.121756 0.124682	0.00771299 -0.121756 0.136947	0.00771299 -0.121756 0.136947	-0.00235252 -0.121756 0.124682	-0.00105944 -0.134885 0.124682	0.00771299 -0.121756 0.136947	-0.00105944 -0.134885 0.124682	0.00881267 -0.132921 0.136947	0.00881267 -0.132921 0.136947	-0.00105944 -0.134885 0.124682	0.00277012 -0.147509 0.124682	0.00881267 -0.132921 0.136947	0.00277012 -0.147509 0.124682	0.0120694 -0.143657 0.136947	0.0120694 -0.143657 0.136947	0.00277012 -0.147509 0.124682	0.00898898 -0.159144 0.124682	0.0120694 -0.143657 0.136947	0.00898898 -0.159144 0.124682	0.0173581 -0.153552 0.136947	0.0173581 -0.153552 0.136947	0.00898898 -0.159144 0.124682	0.0173581 -0.169342 0.124682	0.0173581 -0.153552 0.136947	0.0173581 -0.169342 0.124682	0.0244755 -0.162224 0.136947	0.0244755 -0.162224 0.136947	0.0173581 -0.169342 0.124682	0.027556 -0.177711 0.124682	0.0244755 -0.162224 0.136947	0.027556 -0.177711 0.124682	0.0331481 -0.169342 0.136947	0.0331481 -0.169342 0.136947	0.027556 -0.177711 0.124682	0.0391907 -0.18393 0.124682	0.0331481 -0.169342 0.136947	0.0391907 -0.18393 0.124682	0.0430426 -0.17463 0.136947	0.0430426 -0.17463 0.136947	0.0391907 -0.18393 0.124682	0.051815 -0.187759 0.124682	0.0430426 -0.17463 0.136947	0.051815 -0.187759 0.124682	0.0537787 -0.177887 0.136947	0.0537787 -0.177887 0.136947	0.051815 -0.187759 0.124682	0.0649439 -0.189052 0.124682	0.0537787 -0.177887 0.136947	0.0649439 -0.189052 0.124682	0.0649439 -0.178987 0.136947	0.0649439 -0.178987 0.136947	0.0649439 -0.189052 0.124682	0.0780727 -0.187759 0.124682	0.0649439 -0.178987 0.136947	0.0780727 -0.187759 0.124682	0.0761091 -0.177887 0.136947	0.0761091 -0.177887 0.136947	0.0780727 -0.187759 0.124682	0.0906971 -0.18393 0.124682	0.0761091 -0.177887 0.136947	0.0906971 -0.18393 0.124682	0.0868452 -0.17463 0.136947	0.0868452 -0.17463 0.136947	0.0906971 -0.18393 0.124682	0.102332 -0.177711 0.124682	0.0868452 -0.17463 0.136947	0.102332 -0.177711 0.124682	0.0967396 -0.169342 0.136947	0.0967396 -0.169342 0.136947	0.102332 -0.177711 0.124682	0.11253 -0.169342 0.124682	0.0967396 -0.169342 0.136947	0.11253 -0.169342 0.124682	0.105412 -0.162224 0.136947	0.105412 -0.162224 0.136947	0.11253 -0.169342 0.124682	0.120899 -0.159144 0.124682	0.105412 -0.162224 0.136947	0.120899 -0.159144 0.124682	0.11253 -0.153552 0.136947	0.11253 -0.153552 0.136947	0.120899 -0.159144 0.124682	0.127118 -0.147509 0.124682	0.11253 -0.153552 0.136947	0.127118 -0.147509 0.124682	0.117818 -0.143657 0.136947	0.117818 -0.143657 0.136947	0.127118 -0.147509 0.124682	0.130947 -0.134885 0.124682	0.117818 -0.143657 0.136947	0.130947 -0.134885 0.124682	0.121075 -0.132921 0.136947	0.121075 -0.132921 0.136947	0.130947 -0.134885 0.124682	0.13224 -0.121756 0.124682	0.121075 -0.132921 0.136947	0.13224 -0.121756 0.124682	0.122175 -0.121756 0.136947	0.122175 -0.121756 0.136947	0.13224 -0.121756 0.124682	0.130947 -0.108627 0.124682	0.122175 -0.121756 0.136947	0.130947 -0.108627 0.124682	0.121075 -0.110591 0.136947	0.121075 -0.110591 0.136947	0.130947 -0.108627 0.124682	0.127117 -0.0960028 0.124682	0.121075 -0.110591 0.136947	0.127117 -0.0960028 0.124682	0.117818 -0.0998547 0.136947	0.117818 -0.0998547 0.136947	0.127117 -0.0960028 0.124682	0.120899 -0.0843681 0.124682	0.117818 -0.0998547 0.136947	0.120899 -0.0843681 0.124682	0.112529 -0.0899602 0.136947	0.112529 -0.0899602 0.136947	0.120899 -0.0843681 0.124682	0.112529 -0.0741703 0.124682	0.112529 -0.0899602 0.136947	0.112529 -0.0741703 0.124682	0.105412 -0.0812877 0.136947	0.105412 -0.0812877 0.136947	0.112529 -0.0741703 0.124682	0.102332 -0.0658011 0.124682	0.105412 -0.0812877 0.136947	0.102332 -0.0658011 0.124682	0.0967395 -0.0741703 0.136947	0.0967395 -0.0741703 0.136947	0.102332 -0.0658011 0.124682	0.0906969 -0.0595823 0.124682	0.0967395 -0.0741703 0.136947	0.0906969 -0.0595823 0.124682	0.086845 -0.0688816 0.136947	0.086845 -0.0688816 0.136947	0.0906969 -0.0595823 0.124682	0.0780725 -0.0557528 0.124682	0.086845 -0.0688816 0.136947	0.0780725 -0.0557528 0.124682	0.0761089 -0.0656249 0.136947	0.0761089 -0.0656249 0.136947	0.0780725 -0.0557528 0.124682	0.0649438 -0.0544597 0.124682	0.0761089 -0.0656249 0.136947	0.0649438 -0.0544597 0.124682	0.0649438 -0.0645252 0.136947	0.0649438 -0.0544597 0.124682	0.0649438 -0.0469804 0.110689	0.0503558 -0.0484172 0.110689	0.0649438 -0.0544597 0.124682	0.0503558 -0.0484172 0.110689	0.051815 -0.0557528 0.124682	0.051815 -0.0557528 0.124682	0.0503558 -0.0484172 0.110689	0.0363284 -0.0526723 0.110689	0.051815 -0.0557528 0.124682	0.0363284 -0.0526723 0.110689	0.0391906 -0.0595824 0.124682	0.0391906 -0.0595824 0.124682	0.0363284 -0.0526723 0.110689	0.0234007 -0.0595824 0.110689	0.0391906 -0.0595824 0.124682	0.0234007 -0.0595824 0.110689	0.027556 -0.0658012 0.124682	0.027556 -0.0658012 0.124682	0.0234007 -0.0595824 0.110689	0.0120694 -0.0688817 0.110689	0.027556 -0.0658012 0.124682	0.0120694 -0.0688817 0.110689	0.0173581 -0.0741704 0.124682	0.0173581 -0.0741704 0.124682	0.0120694 -0.0688817 0.110689	0.00277011 -0.0802129 0.110689	0.0173581 -0.0741704 0.124682	0.00277011 -0.0802129 0.110689	0.00898895 -0.0843682 0.124682	0.00898895 -0.0843682 0.124682	0.00277011 -0.0802129 0.110689	-0.00413991 -0.0931407 0.110689	0.00898895 -0.0843682 0.124682	-0.00413991 -0.0931407 0.110689	0.0027701 -0.0960029 0.124682	0.0027701 -0.0960029 0.124682	-0.00413991 -0.0931407 0.110689	-0.00839507 -0.107168 0.110689	0.0027701 -0.0960029 0.124682	-0.00839507 -0.107168 0.110689	-0.00105944 -0.108627 0.124682	-0.00105944 -0.108627 0.124682	-0.00839507 -0.107168 0.110689	-0.00983186 -0.121756 0.110689	-0.00105944 -0.108627 0.124682	-0.00983186 -0.121756 0.110689	-0.00235252 -0.121756 0.124682	-0.00235252 -0.121756 0.124682	-0.00983186 -0.121756 0.110689	-0.00839506 -0.136344 0.110689	-0.00235252 -0.121756 0.124682	-0.00839506 -0.136344 0.110689	-0.00105944 -0.134885 0.124682	-0.00105944 -0.134885 0.124682	-0.00839506 -0.136344 0.110689	-0.00413989 -0.150372 0.110689	-0.00105944 -0.134885 0.124682	-0.00413989 -0.150372 0.110689	0.00277012 -0.147509 0.124682	0.00277012 -0.147509 0.124682	-0.00413989 -0.150372 0.110689	0.00277014 -0.163299 0.110689	0.00277012 -0.147509 0.124682	0.00277014 -0.163299 0.110689	0.00898898 -0.159144 0.124682	0.00898898 -0.159144 0.124682	0.00277014 -0.163299 0.110689	0.0120695 -0.17463 0.110689	0.00898898 -0.159144 0.124682	0.0120695 -0.17463 0.110689	0.0173581 -0.169342 0.124682	0.0173581 -0.169342 0.124682	0.0120695 -0.17463 0.110689	0.0234007 -0.18393 0.110689	0.0173581 -0.169342 0.124682	0.0234007 -0.18393 0.110689	0.027556 -0.177711 0.124682	0.027556 -0.177711 0.124682	0.0234007 -0.18393 0.110689	0.0363285 -0.19084 0.110689	0.027556 -0.177711 0.124682	0.0363285 -0.19084 0.110689	0.0391907 -0.18393 0.124682	0.0391907 -0.18393 0.124682	0.0363285 -0.19084 0.110689	0.0503559 -0.195095 0.110689	0.0391907 -0.18393 0.124682	0.0503559 -0.195095 0.110689	0.051815 -0.187759 0.124682	0.051815 -0.187759 0.124682	0.0503559 -0.195095 0.110689	0.0649439 -0.196532 0.110689	0.051815 -0.187759 0.124682	0.0649439 -0.196532 0.110689	0.0649439 -0.189052 0.124682	0.0649439 -0.189052 0.124682	0.0649439 -0.196532 0.110689	0.0795319 -0.195095 0.110689	0.0649439 -0.189052 0.124682	0.0795319 -0.195095 0.110689	0.0780727 -0.187759 0.124682	0.0780727 -0.187759 0.124682	0.0795319 -0.195095 0.110689	0.0935593 -0.19084 0.110689	0.0780727 -0.187759 0.124682	0.0935593 -0.19084 0.110689	0.0906971 -0.18393 0.124682	0.0906971 -0.18393 0.124682	0.0935593 -0.19084 0.110689	0.106487 -0.18393 0.110689	0.0906971 -0.18393 0.124682	0.106487 -0.18393 0.110689	0.102332 -0.177711 0.124682	0.102332 -0.177711 0.124682	0.106487 -0.18393 0.110689	0.117818 -0.17463 0.110689	0.102332 -0.177711 0.124682	0.117818 -0.17463 0.110689	0.11253 -0.169342 0.124682	0.11253 -0.169342 0.124682	0.117818 -0.17463 0.110689	0.127118 -0.163299 0.110689	0.11253 -0.169342 0.124682	0.127118 -0.163299 0.110689	0.120899 -0.159144 0.124682	0.120899 -0.159144 0.124682	0.127118 -0.163299 0.110689	0.134028 -0.150371 0.110689	0.120899 -0.159144 0.124682	0.134028 -0.150371 0.110689	0.127118 -0.147509 0.124682	0.127118 -0.147509 0.124682	0.134028 -0.150371 0.110689	0.138283 -0.136344 0.110689	0.127118 -0.147509 0.124682	0.138283 -0.136344 0.110689	0.130947 -0.134885 0.124682	0.130947 -0.134885 0.124682	0.138283 -0.136344 0.110689	0.139719 -0.121756 0.110689	0.130947 -0.134885 0.124682	0.139719 -0.121756 0.110689	0.13224 -0.121756 0.124682	0.13224 -0.121756 0.124682	0.139719 -0.121756 0.110689	0.138283 -0.107168 0.110689	0.13224 -0.121756 0.124682	0.138283 -0.107168 0.110689	0.130947 -0.108627 0.124682	0.130947 -0.108627 0.124682	0.138283 -0.107168 0.110689	0.134027 -0.0931405 0.110689	0.130947 -0.108627 0.124682	0.134027 -0.0931405 0.110689	0.127117 -0.0960028 0.124682	0.127117 -0.0960028 0.124682	0.134027 -0.0931405 0.110689	0.127117 -0.0802128 0.110689	0.127117 -0.0960028 0.124682	0.127117 -0.0802128 0.110689	0.120899 -0.0843681 0.124682	0.120899 -0.0843681 0.124682	0.127117 -0.0802128 0.110689	0.117818 -0.0688816 0.110689	0.120899 -0.0843681 0.124682	0.117818 -0.0688816 0.110689	0.112529 -0.0741703 0.124682	0.112529 -0.0741703 0.124682	0.117818 -0.0688816 0.110689	0.106487 -0.0595823 0.110689	0.112529 -0.0741703 0.124682	0.106487 -0.0595823 0.110689	0.102332 -0.0658011 0.124682	0.102332 -0.0658011 0.124682	0.106487 -0.0595823 0.110689	0.0935591 -0.0526723 0.110689	0.102332 -0.0658011 0.124682	0.0935591 -0.0526723 0.110689	0.0906969 -0.0595823 0.124682	0.0906969 -0.0595823 0.124682	0.0935591 -0.0526723 0.110689	0.0795317 -0.0484172 0.110689	0.0906969 -0.0595823 0.124682	0.0795317 -0.0484172 0.110689	0.0780725 -0.0557528 0.124682	0.0780725 -0.0557528 0.124682	0.0795317 -0.0484172 0.110689	0.0649438 -0.0469804 0.110689	0.0780725 -0.0557528 0.124682	0.0649438 -0.0469804 0.110689	0.0649438 -0.0544597 0.124682	0.0649438 -0.0469804 0.110689	0.0649438 -0.0423746 0.0955061	0.0494573 -0.0438999 0.0955061	0.0649438 -0.0469804 0.110689	0.0494573 -0.0438999 0.0955061	0.0503558 -0.0484172 0.110689	0.0503558 -0.0484172 0.110689	0.0494573 -0.0438999 0.0955061	0.0345659 -0.0484172 0.0955061	0.0503558 -0.0484172 0.110689	0.0345659 -0.0484172 0.0955061	0.0363284 -0.0526723 0.110689	0.0363284 -0.0526723 0.110689	0.0345659 -0.0484172 0.0955061	0.0208419 -0.0557528 0.0955061	0.0363284 -0.0526723 0.110689	0.0208419 -0.0557528 0.0955061	0.0234007 -0.0595824 0.110689	0.0234007 -0.0595824 0.110689	0.0208419 -0.0557528 0.0955061	0.00881267 -0.0656249 0.0955061	0.0234007 -0.0595824 0.110689	0.00881267 -0.0656249 0.0955061	0.0120694 -0.0688817 0.110689	0.0120694 -0.0688817 0.110689	0.00881267 -0.0656249 0.0955061	-0.00105944 -0.0776541 0.0955061	0.0120694 -0.0688817 0.110689	-0.00105944 -0.0776541 0.0955061	0.00277011 -0.0802129 0.110689	0.00277011 -0.0802129 0.110689	-0.00105944 -0.0776541 0.0955061	-0.00839507 -0.0913781 0.0955061	0.00277011 -0.0802129 0.110689	-0.00839507 -0.0913781 0.0955061	-0.00413991 -0.0931407 0.110689	-0.00413991 -0.0931407 0.110689	-0.00839507 -0.0913781 0.0955061	-0.0129123 -0.10627 0.0955061	-0.00413991 -0.0931407 0.110689	-0.0129123 -0.10627 0.0955061	-0.00839507 -0.107168 0.110689	-0.00839507 -0.107168 0.110689	-0.0129123 -0.10627 0.0955061	-0.0144376 -0.121756 0.0955061	-0.00839507 -0.107168 0.110689	-0.0144376 -0.121756 0.0955061	-0.00983186 -0.121756 0.110689	-0.00983186 -0.121756 0.110689	-0.0144376 -0.121756 0.0955061	-0.0129123 -0.137243 0.0955061	-0.00983186 -0.121756 0.110689	-0.0129123 -0.137243 0.0955061	-0.00839506 -0.136344 0.110689	-0.00839506 -0.136344 0.110689	-0.0129123 -0.137243 0.0955061	-0.00839505 -0.152134 0.0955061	-0.00839506 -0.136344 0.110689	-0.00839505 -0.152134 0.0955061	-0.00413989 -0.150372 0.110689	-0.00413989 -0.150372 0.110689	-0.00839505 -0.152134 0.0955061	-0.00105941 -0.165858 0.0955061	-0.00413989 -0.150372 0.110689	-0.00105941 -0.165858 0.0955061	0.00277014 -0.163299 0.110689	0.00277014 -0.163299 0.110689	-0.00105941 -0.165858 0.0955061	0.0088127 -0.177887 0.0955061	0.00277014 -0.163299 0.110689	0.0088127 -0.177887 0.0955061	0.0120695 -0.17463 0.110689	0.0120695 -0.17463 0.110689	0.0088127 -0.177887 0.0955061	0.0208419 -0.187759 0.0955061	0.0120695 -0.17463 0.110689	0.0208419 -0.187759 0.0955061	0.0234007 -0.18393 0.110689	0.0234007 -0.18393 0.110689	0.0208419 -0.187759 0.0955061	0.0345659 -0.195095 0.0955061	0.0234007 -0.18393 0.110689	0.0345659 -0.195095 0.0955061	0.0363285 -0.19084 0.110689	0.0363285 -0.19084 0.110689	0.0345659 -0.195095 0.0955061	0.0494573 -0.199612 0.0955061	0.0363285 -0.19084 0.110689	0.0494573 -0.199612 0.0955061	0.0503559 -0.195095 0.110689	0.0503559 -0.195095 0.110689	0.0494573 -0.199612 0.0955061	0.0649439 -0.201137 0.0955061	0.0503559 -0.195095 0.110689	0.0649439 -0.201137 0.0955061	0.0649439 -0.196532 0.110689	0.0649439 -0.196532 0.110689	0.0649439 -0.201137 0.0955061	0.0804304 -0.199612 0.0955061	0.0649439 -0.196532 0.110689	0.0804304 -0.199612 0.0955061	0.0795319 -0.195095 0.110689	0.0795319 -0.195095 0.110689	0.0804304 -0.199612 0.0955061	0.0953218 -0.195095 0.0955061	0.0795319 -0.195095 0.110689	0.0953218 -0.195095 0.0955061	0.0935593 -0.19084 0.110689	0.0935593 -0.19084 0.110689	0.0953218 -0.195095 0.0955061	0.109046 -0.187759 0.0955061	0.0935593 -0.19084 0.110689	0.109046 -0.187759 0.0955061	0.106487 -0.18393 0.110689	0.106487 -0.18393 0.110689	0.109046 -0.187759 0.0955061	0.121075 -0.177887 0.0955061	0.106487 -0.18393 0.110689	0.121075 -0.177887 0.0955061	0.117818 -0.17463 0.110689	0.117818 -0.17463 0.110689	0.121075 -0.177887 0.0955061	0.130947 -0.165858 0.0955061	0.117818 -0.17463 0.110689	0.130947 -0.165858 0.0955061	0.127118 -0.163299 0.110689	0.127118 -0.163299 0.110689	0.130947 -0.165858 0.0955061	0.138283 -0.152134 0.0955061	0.127118 -0.163299 0.110689	0.138283 -0.152134 0.0955061	0.134028 -0.150371 0.110689	0.134028 -0.150371 0.110689	0.138283 -0.152134 0.0955061	0.1428 -0.137242 0.0955061	0.134028 -0.150371 0.110689	0.1428 -0.137242 0.0955061	0.138283 -0.136344 0.110689	0.138283 -0.136344 0.110689	0.1428 -0.137242 0.0955061	0.144325 -0.121756 0.0955061	0.138283 -0.136344 0.110689	0.144325 -0.121756 0.0955061	0.139719 -0.121756 0.110689	0.139719 -0.121756 0.110689	0.144325 -0.121756 0.0955061	0.1428 -0.106269 0.0955061	0.139719 -0.121756 0.110689	0.1428 -0.106269 0.0955061	0.138283 -0.107168 0.110689	0.138283 -0.107168 0.110689	0.1428 -0.106269 0.0955061	0.138283 -0.091378 0.0955061	0.138283 -0.107168 0.110689	0.138283 -0.091378 0.0955061	0.134027 -0.0931405 0.110689	0.134027 -0.0931405 0.110689	0.138283 -0.091378 0.0955061	0.130947 -0.077654 0.0955061	0.134027 -0.0931405 0.110689	0.130947 -0.077654 0.0955061	0.127117 -0.0802128 0.110689	0.127117 -0.0802128 0.110689	0.130947 -0.077654 0.0955061	0.121075 -0.0656248 0.0955061	0.127117 -0.0802128 0.110689	0.121075 -0.0656248 0.0955061	0.117818 -0.0688816 0.110689	0.117818 -0.0688816 0.110689	0.121075 -0.0656248 0.0955061	0.109046 -0.0557527 0.0955061	0.117818 -0.0688816 0.110689	0.109046 -0.0557527 0.0955061	0.106487 -0.0595823 0.110689	0.106487 -0.0595823 0.110689	0.109046 -0.0557527 0.0955061	0.0953216 -0.0484171 0.0955061	0.106487 -0.0595823 0.110689	0.0953216 -0.0484171 0.0955061	0.0935591 -0.0526723 0.110689	0.0935591 -0.0526723 0.110689	0.0953216 -0.0484171 0.0955061	0.0804302 -0.0438999 0.0955061	0.0935591 -0.0526723 0.110689	0.0804302 -0.0438999 0.0955061	0.0795317 -0.0484172 0.110689	0.0795317 -0.0484172 0.110689	0.0804302 -0.0438999 0.0955061	0.0649438 -0.0423746 0.0955061	0.0795317 -0.0484172 0.110689	0.0649438 -0.0423746 0.0955061	0.0649438 -0.0469804 0.110689	0.0649438 -0.0423746 0.0955061	0.0649438 -0.0408195 0.0797162	0.0491539 -0.0423746 0.0797162	0.0649438 -0.0423746 0.0955061	0.0491539 -0.0423746 0.0797162	0.0494573 -0.0438999 0.0955061	0.0494573 -0.0438999 0.0955061	0.0491539 -0.0423746 0.0797162	0.0339707 -0.0469804 0.0797162	0.0494573 -0.0438999 0.0955061	0.0339707 -0.0469804 0.0797162	0.0345659 -0.0484172 0.0955061	0.0345659 -0.0484172 0.0955061	0.0339707 -0.0469804 0.0797162	0.0199778 -0.0544597 0.0797162	0.0345659 -0.0484172 0.0955061	0.0199778 -0.0544597 0.0797162	0.0208419 -0.0557528 0.0955061	0.0208419 -0.0557528 0.0955061	0.0199778 -0.0544597 0.0797162	0.00771299 -0.0645252 0.0797162	0.0208419 -0.0557528 0.0955061	0.00771299 -0.0645252 0.0797162	0.00881267 -0.0656249 0.0955061	0.00881267 -0.0656249 0.0955061	0.00771299 -0.0645252 0.0797162	-0.00235252 -0.0767901 0.0797162	0.00881267 -0.0656249 0.0955061	-0.00235252 -0.0767901 0.0797162	-0.00105944 -0.0776541 0.0955061	-0.00105944 -0.0776541 0.0955061	-0.00235252 -0.0767901 0.0797162	-0.00983187 -0.090783 0.0797162	-0.00105944 -0.0776541 0.0955061	-0.00983187 -0.090783 0.0797162	-0.00839507 -0.0913781 0.0955061	-0.00839507 -0.0913781 0.0955061	-0.00983187 -0.090783 0.0797162	-0.0144376 -0.105966 0.0797162	-0.00839507 -0.0913781 0.0955061	-0.0144376 -0.105966 0.0797162	-0.0129123 -0.10627 0.0955061	-0.0129123 -0.10627 0.0955061	-0.0144376 -0.105966 0.0797162	-0.0159928 -0.121756 0.0797162	-0.0129123 -0.10627 0.0955061	-0.0159928 -0.121756 0.0797162	-0.0144376 -0.121756 0.0955061	-0.0144376 -0.121756 0.0955061	-0.0159928 -0.121756 0.0797162	-0.0144376 -0.137546 0.0797162	-0.0144376 -0.121756 0.0955061	-0.0144376 -0.137546 0.0797162	-0.0129123 -0.137243 0.0955061	-0.0129123 -0.137243 0.0955061	-0.0144376 -0.137546 0.0797162	-0.00983185 -0.152729 0.0797162	-0.0129123 -0.137243 0.0955061	-0.00983185 -0.152729 0.0797162	-0.00839505 -0.152134 0.0955061	-0.00839505 -0.152134 0.0955061	-0.00983185 -0.152729 0.0797162	-0.00235249 -0.166722 0.0797162	-0.00839505 -0.152134 0.0955061	-0.00235249 -0.166722 0.0797162	-0.00105941 -0.165858 0.0955061	-0.00105941 -0.165858 0.0955061	-0.00235249 -0.166722 0.0797162	0.00771302 -0.178987 0.0797162	-0.00105941 -0.165858 0.0955061	0.00771302 -0.178987 0.0797162	0.0088127 -0.177887 0.0955061	0.0088127 -0.177887 0.0955061	0.00771302 -0.178987 0.0797162	0.0199779 -0.189052 0.0797162	0.0088127 -0.177887 0.0955061	0.0199779 -0.189052 0.0797162	0.0208419 -0.187759 0.0955061	0.0208419 -0.187759 0.0955061	0.0199779 -0.189052 0.0797162	0.0339708 -0.196532 0.0797162	0.0208419 -0.187759 0.0955061	0.0339708 -0.196532 0.0797162	0.0345659 -0.195095 0.0955061	0.0345659 -0.195095 0.0955061	0.0339708 -0.196532 0.0797162	0.0491539 -0.201138 0.0797162	0.0345659 -0.195095 0.0955061	0.0491539 -0.201138 0.0797162	0.0494573 -0.199612 0.0955061	0.0494573 -0.199612 0.0955061	0.0491539 -0.201138 0.0797162	0.0649439 -0.202693 0.0797162	0.0494573 -0.199612 0.0955061	0.0649439 -0.202693 0.0797162	0.0649439 -0.201137 0.0955061	0.0649439 -0.201137 0.0955061	0.0649439 -0.202693 0.0797162	0.0807338 -0.201137 0.0797162	0.0649439 -0.201137 0.0955061	0.0807338 -0.201137 0.0797162	0.0804304 -0.199612 0.0955061	0.0804304 -0.199612 0.0955061	0.0807338 -0.201137 0.0797162	0.095917 -0.196532 0.0797162	0.0804304 -0.199612 0.0955061	0.095917 -0.196532 0.0797162	0.0953218 -0.195095 0.0955061	0.0953218 -0.195095 0.0955061	0.095917 -0.196532 0.0797162	0.10991 -0.189052 0.0797162	0.0953218 -0.195095 0.0955061	0.10991 -0.189052 0.0797162	0.109046 -0.187759 0.0955061	0.109046 -0.187759 0.0955061	0.10991 -0.189052 0.0797162	0.122175 -0.178987 0.0797162	0.109046 -0.187759 0.0955061	0.122175 -0.178987 0.0797162	0.121075 -0.177887 0.0955061	0.121075 -0.177887 0.0955061	0.122175 -0.178987 0.0797162	0.13224 -0.166722 0.0797162	0.121075 -0.177887 0.0955061	0.13224 -0.166722 0.0797162	0.130947 -0.165858 0.0955061	0.130947 -0.165858 0.0955061	0.13224 -0.166722 0.0797162	0.13972 -0.152729 0.0797162	0.130947 -0.165858 0.0955061	0.13972 -0.152729 0.0797162	0.138283 -0.152134 0.0955061	0.138283 -0.152134 0.0955061	0.13972 -0.152729 0.0797162	0.144325 -0.137546 0.0797162	0.138283 -0.152134 0.0955061	0.144325 -0.137546 0.0797162	0.1428 -0.137242 0.0955061	0.1428 -0.137242 0.0955061	0.144325 -0.137546 0.0797162	0.14588 -0.121756 0.0797162	0.1428 -0.137242 0.0955061	0.14588 -0.121756 0.0797162	0.144325 -0.121756 0.0955061	0.144325 -0.121756 0.0955061	0.14588 -0.121756 0.0797162	0.144325 -0.105966 0.0797162	0.144325 -0.121756 0.0955061	0.144325 -0.105966 0.0797162	0.1428 -0.106269 0.0955061	0.1428 -0.106269 0.0955061	0.144325 -0.105966 0.0797162	0.139719 -0.0907828 0.0797162	0.1428 -0.106269 0.0955061	0.139719 -0.0907828 0.0797162	0.138283 -0.091378 0.0955061	0.138283 -0.091378 0.0955061	0.139719 -0.0907828 0.0797162	0.13224 -0.07679 0.0797162	0.138283 -0.091378 0.0955061	0.13224 -0.07679 0.0797162	0.130947 -0.077654 0.0955061	0.130947 -0.077654 0.0955061	0.13224 -0.07679 0.0797162	0.122175 -0.0645251 0.0797162	0.130947 -0.077654 0.0955061	0.122175 -0.0645251 0.0797162	0.121075 -0.0656248 0.0955061	0.121075 -0.0656248 0.0955061	0.122175 -0.0645251 0.0797162	0.10991 -0.0544596 0.0797162	0.121075 -0.0656248 0.0955061	0.10991 -0.0544596 0.0797162	0.109046 -0.0557527 0.0955061	0.109046 -0.0557527 0.0955061	0.10991 -0.0544596 0.0797162	0.0959168 -0.0469803 0.0797162	0.109046 -0.0557527 0.0955061	0.0959168 -0.0469803 0.0797162	0.0953216 -0.0484171 0.0955061	0.0953216 -0.0484171 0.0955061	0.0959168 -0.0469803 0.0797162	0.0807336 -0.0423746 0.0797162	0.0953216 -0.0484171 0.0955061	0.0807336 -0.0423746 0.0797162	0.0804302 -0.0438999 0.0955061	0.0804302 -0.0438999 0.0955061	0.0807336 -0.0423746 0.0797162	0.0649438 -0.0408195 0.0797162	0.0804302 -0.0438999 0.0955061	0.0649438 -0.0408195 0.0797162	0.0649438 -0.0423746 0.0955061	0.0649438 -0.0408195 0.0797162	0.0649438 -0.0423746 0.0639262	0.0494573 -0.0438999 0.0639262	0.0649438 -0.0408195 0.0797162	0.0494573 -0.0438999 0.0639262	0.0491539 -0.0423746 0.0797162	0.0491539 -0.0423746 0.0797162	0.0494573 -0.0438999 0.0639262	0.0345659 -0.0484172 0.0639262	0.0491539 -0.0423746 0.0797162	0.0345659 -0.0484172 0.0639262	0.0339707 -0.0469804 0.0797162	0.0339707 -0.0469804 0.0797162	0.0345659 -0.0484172 0.0639262	0.0208419 -0.0557528 0.0639262	0.0339707 -0.0469804 0.0797162	0.0208419 -0.0557528 0.0639262	0.0199778 -0.0544597 0.0797162	0.0199778 -0.0544597 0.0797162	0.0208419 -0.0557528 0.0639262	0.00881266 -0.0656249 0.0639262	0.0199778 -0.0544597 0.0797162	0.00881266 -0.0656249 0.0639262	0.00771299 -0.0645252 0.0797162	0.00771299 -0.0645252 0.0797162	0.00881266 -0.0656249 0.0639262	-0.00105944 -0.0776541 0.0639262	0.00771299 -0.0645252 0.0797162	-0.00105944 -0.0776541 0.0639262	-0.00235252 -0.0767901 0.0797162	-0.00235252 -0.0767901 0.0797162	-0.00105944 -0.0776541 0.0639262	-0.00839508 -0.0913781 0.0639262	-0.00235252 -0.0767901 0.0797162	-0.00839508 -0.0913781 0.0639262	-0.00983187 -0.090783 0.0797162	-0.00983187 -0.090783 0.0797162	-0.00839508 -0.0913781 0.0639262	-0.0129123 -0.10627 0.0639262	-0.00983187 -0.090783 0.0797162	-0.0129123 -0.10627 0.0639262	-0.0144376 -0.105966 0.0797162	-0.0144376 -0.105966 0.0797162	-0.0129123 -0.10627 0.0639262	-0.0144376 -0.121756 0.0639262	-0.0144376 -0.105966 0.0797162	-0.0144376 -0.121756 0.0639262	-0.0159928 -0.121756 0.0797162	-0.0159928 -0.121756 0.0797162	-0.0144376 -0.121756 0.0639262	-0.0129123 -0.137243 0.0639262	-0.0159928 -0.121756 0.0797162	-0.0129123 -0.137243 0.0639262	-0.0144376 -0.137546 0.0797162	-0.0144376 -0.137546 0.0797162	-0.0129123 -0.137243 0.0639262	-0.00839505 -0.152134 0.0639262	-0.0144376 -0.137546 0.0797162	-0.00839505 -0.152134 0.0639262	-0.00983185 -0.152729 0.0797162	-0.00983185 -0.152729 0.0797162	-0.00839505 -0.152134 0.0639262	-0.00105941 -0.165858 0.0639262	-0.00983185 -0.152729 0.0797162	-0.00105941 -0.165858 0.0639262	-0.00235249 -0.166722 0.0797162	-0.00235249 -0.166722 0.0797162	-0.00105941 -0.165858 0.0639262	0.00881269 -0.177887 0.0639262	-0.00235249 -0.166722 0.0797162	0.00881269 -0.177887 0.0639262	0.00771302 -0.178987 0.0797162	0.00771302 -0.178987 0.0797162	0.00881269 -0.177887 0.0639262	0.0208419 -0.187759 0.0639262	0.00771302 -0.178987 0.0797162	0.0208419 -0.187759 0.0639262	0.0199779 -0.189052 0.0797162	0.0199779 -0.189052 0.0797162	0.0208419 -0.187759 0.0639262	0.0345659 -0.195095 0.0639262	0.0199779 -0.189052 0.0797162	0.0345659 -0.195095 0.0639262	0.0339708 -0.196532 0.0797162	0.0339708 -0.196532 0.0797162	0.0345659 -0.195095 0.0639262	0.0494573 -0.199612 0.0639262	0.0339708 -0.196532 0.0797162	0.0494573 -0.199612 0.0639262	0.0491539 -0.201138 0.0797162	0.0491539 -0.201138 0.0797162	0.0494573 -0.199612 0.0639262	0.0649439 -0.201137 0.0639262	0.0491539 -0.201138 0.0797162	0.0649439 -0.201137 0.0639262	0.0649439 -0.202693 0.0797162	0.0649439 -0.202693 0.0797162	0.0649439 -0.201137 0.0639262	0.0804304 -0.199612 0.0639262	0.0649439 -0.202693 0.0797162	0.0804304 -0.199612 0.0639262	0.0807338 -0.201137 0.0797162	0.0807338 -0.201137 0.0797162	0.0804304 -0.199612 0.0639262	0.0953219 -0.195095 0.0639262	0.0807338 -0.201137 0.0797162	0.0953219 -0.195095 0.0639262	0.095917 -0.196532 0.0797162	0.095917 -0.196532 0.0797162	0.0953219 -0.195095 0.0639262	0.109046 -0.187759 0.0639262	0.095917 -0.196532 0.0797162	0.109046 -0.187759 0.0639262	0.10991 -0.189052 0.0797162	0.10991 -0.189052 0.0797162	0.109046 -0.187759 0.0639262	0.121075 -0.177887 0.0639262	0.10991 -0.189052 0.0797162	0.121075 -0.177887 0.0639262	0.122175 -0.178987 0.0797162	0.122175 -0.178987 0.0797162	0.121075 -0.177887 0.0639262	0.130947 -0.165858 0.0639262	0.122175 -0.178987 0.0797162	0.130947 -0.165858 0.0639262	0.13224 -0.166722 0.0797162	0.13224 -0.166722 0.0797162	0.130947 -0.165858 0.0639262	0.138283 -0.152134 0.0639262	0.13224 -0.166722 0.0797162	0.138283 -0.152134 0.0639262	0.13972 -0.152729 0.0797162	0.13972 -0.152729 0.0797162	0.138283 -0.152134 0.0639262	0.1428 -0.137242 0.0639262	0.13972 -0.152729 0.0797162	0.1428 -0.137242 0.0639262	0.144325 -0.137546 0.0797162	0.144325 -0.137546 0.0797162	0.1428 -0.137242 0.0639262	0.144325 -0.121756 0.0639262	0.144325 -0.137546 0.0797162	0.144325 -0.121756 0.0639262	0.14588 -0.121756 0.0797162	0.14588 -0.121756 0.0797162	0.144325 -0.121756 0.0639262	0.1428 -0.106269 0.0639262	0.14588 -0.121756 0.0797162	0.1428 -0.106269 0.0639262	0.144325 -0.105966 0.0797162	0.144325 -0.105966 0.0797162	0.1428 -0.106269 0.0639262	0.138283 -0.091378 0.0639262	0.144325 -0.105966 0.0797162	0.138283 -0.091378 0.0639262	0.139719 -0.0907828 0.0797162	0.139719 -0.0907828 0.0797162	0.138283 -0.091378 0.0639262	0.130947 -0.077654 0.0639262	0.139719 -0.0907828 0.0797162	0.130947 -0.077654 0.0639262	0.13224 -0.07679 0.0797162	0.13224 -0.07679 0.0797162	0.130947 -0.077654 0.0639262	0.121075 -0.0656248 0.0639262	0.13224 -0.07679 0.0797162	0.121075 -0.0656248 0.0639262	0.122175 -0.0645251 0.0797162	0.122175 -0.0645251 0.0797162	0.121075 -0.0656248 0.0639262	0.109046 -0.0557527 0.0639262	0.122175 -0.0645251 0.0797162	0.109046 -0.0557527 0.0639262	0.10991 -0.0544596 0.0797162	0.10991 -0.0544596 0.0797162	0.109046 -0.0557527 0.0639262	0.0953216 -0.0484171 0.0639262	0.10991 -0.0544596 0.0797162	0.0953216 -0.0484171 0.0639262	0.0959168 -0.0469803 0.0797162	0.0959168 -0.0469803 0.0797162	0.0953216 -0.0484171 0.0639262	0.0804302 -0.0438999 0.0639262	0.0959168 -0.0469803 0.0797162	0.0804302 -0.0438999 0.0639262	0.0807336 -0.0423746 0.0797162	0.0807336 -0.0423746 0.0797162	0.0804302 -0.0438999 0.0639262	0.0649438 -0.0423746 0.0639262	0.0807336 -0.0423746 0.0797162	0.0649438 -0.0423746 0.0639262	0.0649438 -0.0408195 0.0797162	0.0649438 -0.0423746 0.0639262	0.0649438 -0.0469804 0.0487431	0.0503558 -0.0484172 0.0487431	0.0649438 -0.0423746 0.0639262	0.0503558 -0.0484172 0.0487431	0.0494573 -0.0438999 0.0639262	0.0494573 -0.0438999 0.0639262	0.0503558 -0.0484172 0.0487431	0.0363284 -0.0526723 0.0487431	0.0494573 -0.0438999 0.0639262	0.0363284 -0.0526723 0.0487431	0.0345659 -0.0484172 0.0639262	0.0345659 -0.0484172 0.0639262	0.0363284 -0.0526723 0.0487431	0.0234007 -0.0595824 0.0487431	0.0345659 -0.0484172 0.0639262	0.0234007 -0.0595824 0.0487431	0.0208419 -0.0557528 0.0639262	0.0208419 -0.0557528 0.0639262	0.0234007 -0.0595824 0.0487431	0.0120694 -0.0688817 0.0487431	0.0208419 -0.0557528 0.0639262	0.0120694 -0.0688817 0.0487431	0.00881266 -0.0656249 0.0639262	0.00881266 -0.0656249 0.0639262	0.0120694 -0.0688817 0.0487431	0.00277011 -0.0802129 0.0487431	0.00881266 -0.0656249 0.0639262	0.00277011 -0.0802129 0.0487431	-0.00105944 -0.0776541 0.0639262	-0.00105944 -0.0776541 0.0639262	0.00277011 -0.0802129 0.0487431	-0.00413991 -0.0931407 0.0487431	-0.00105944 -0.0776541 0.0639262	-0.00413991 -0.0931407 0.0487431	-0.00839508 -0.0913781 0.0639262	-0.00839508 -0.0913781 0.0639262	-0.00413991 -0.0931407 0.0487431	-0.00839507 -0.107168 0.0487431	-0.00839508 -0.0913781 0.0639262	-0.00839507 -0.107168 0.0487431	-0.0129123 -0.10627 0.0639262	-0.0129123 -0.10627 0.0639262	-0.00839507 -0.107168 0.0487431	-0.00983187 -0.121756 0.0487431	-0.0129123 -0.10627 0.0639262	-0.00983187 -0.121756 0.0487431	-0.0144376 -0.121756 0.0639262	-0.0144376 -0.121756 0.0639262	-0.00983187 -0.121756 0.0487431	-0.00839507 -0.136344 0.0487431	-0.0144376 -0.121756 0.0639262	-0.00839507 -0.136344 0.0487431	-0.0129123 -0.137243 0.0639262	-0.0129123 -0.137243 0.0639262	-0.00839507 -0.136344 0.0487431	-0.00413989 -0.150372 0.0487431	-0.0129123 -0.137243 0.0639262	-0.00413989 -0.150372 0.0487431	-0.00839505 -0.152134 0.0639262	-0.00839505 -0.152134 0.0639262	-0.00413989 -0.150372 0.0487431	0.00277013 -0.163299 0.0487431	-0.00839505 -0.152134 0.0639262	0.00277013 -0.163299 0.0487431	-0.00105941 -0.165858 0.0639262	-0.00105941 -0.165858 0.0639262	0.00277013 -0.163299 0.0487431	0.0120695 -0.17463 0.0487431	-0.00105941 -0.165858 0.0639262	0.0120695 -0.17463 0.0487431	0.00881269 -0.177887 0.0639262	0.00881269 -0.177887 0.0639262	0.0120695 -0.17463 0.0487431	0.0234007 -0.18393 0.0487431	0.00881269 -0.177887 0.0639262	0.0234007 -0.18393 0.0487431	0.0208419 -0.187759 0.0639262	0.0208419 -0.187759 0.0639262	0.0234007 -0.18393 0.0487431	0.0363285 -0.19084 0.0487431	0.0208419 -0.187759 0.0639262	0.0363285 -0.19084 0.0487431	0.0345659 -0.195095 0.0639262	0.0345659 -0.195095 0.0639262	0.0363285 -0.19084 0.0487431	0.0503559 -0.195095 0.0487431	0.0345659 -0.195095 0.0639262	0.0503559 -0.195095 0.0487431	0.0494573 -0.199612 0.0639262	0.0494573 -0.199612 0.0639262	0.0503559 -0.195095 0.0487431	0.0649439 -0.196532 0.0487431	0.0494573 -0.199612 0.0639262	0.0649439 -0.196532 0.0487431	0.0649439 -0.201137 0.0639262	0.0649439 -0.201137 0.0639262	0.0649439 -0.196532 0.0487431	0.0795319 -0.195095 0.0487431	0.0649439 -0.201137 0.0639262	0.0795319 -0.195095 0.0487431	0.0804304 -0.199612 0.0639262	0.0804304 -0.199612 0.0639262	0.0795319 -0.195095 0.0487431	0.0935593 -0.19084 0.0487431	0.0804304 -0.199612 0.0639262	0.0935593 -0.19084 0.0487431	0.0953219 -0.195095 0.0639262	0.0953219 -0.195095 0.0639262	0.0935593 -0.19084 0.0487431	0.106487 -0.18393 0.0487431	0.0953219 -0.195095 0.0639262	0.106487 -0.18393 0.0487431	0.109046 -0.187759 0.0639262	0.109046 -0.187759 0.0639262	0.106487 -0.18393 0.0487431	0.117818 -0.17463 0.0487431	0.109046 -0.187759 0.0639262	0.117818 -0.17463 0.0487431	0.121075 -0.177887 0.0639262	0.121075 -0.177887 0.0639262	0.117818 -0.17463 0.0487431	0.127118 -0.163299 0.0487431	0.121075 -0.177887 0.0639262	0.127118 -0.163299 0.0487431	0.130947 -0.165858 0.0639262	0.130947 -0.165858 0.0639262	0.127118 -0.163299 0.0487431	0.134028 -0.150371 0.0487431	0.130947 -0.165858 0.0639262	0.134028 -0.150371 0.0487431	0.138283 -0.152134 0.0639262	0.138283 -0.152134 0.0639262	0.134028 -0.150371 0.0487431	0.138283 -0.136344 0.0487431	0.138283 -0.152134 0.0639262	0.138283 -0.136344 0.0487431	0.1428 -0.137242 0.0639262	0.1428 -0.137242 0.0639262	0.138283 -0.136344 0.0487431	0.139719 -0.121756 0.0487431	0.1428 -0.137242 0.0639262	0.139719 -0.121756 0.0487431	0.144325 -0.121756 0.0639262	0.144325 -0.121756 0.0639262	0.139719 -0.121756 0.0487431	0.138283 -0.107168 0.0487431	0.144325 -0.121756 0.0639262	0.138283 -0.107168 0.0487431	0.1428 -0.106269 0.0639262	0.1428 -0.106269 0.0639262	0.138283 -0.107168 0.0487431	0.134027 -0.0931405 0.0487431	0.1428 -0.106269 0.0639262	0.134027 -0.0931405 0.0487431	0.138283 -0.091378 0.0639262	0.138283 -0.091378 0.0639262	0.134027 -0.0931405 0.0487431	0.127117 -0.0802128 0.0487431	0.138283 -0.091378 0.0639262	0.127117 -0.0802128 0.0487431	0.130947 -0.077654 0.0639262	0.130947 -0.077654 0.0639262	0.127117 -0.0802128 0.0487431	0.117818 -0.0688816 0.0487431	0.130947 -0.077654 0.0639262	0.117818 -0.0688816 0.0487431	0.121075 -0.0656248 0.0639262	0.121075 -0.0656248 0.0639262	0.117818 -0.0688816 0.0487431	0.106487 -0.0595823 0.0487431	0.121075 -0.0656248 0.0639262	0.106487 -0.0595823 0.0487431	0.109046 -0.0557527 0.0639262	0.109046 -0.0557527 0.0639262	0.106487 -0.0595823 0.0487431	0.0935591 -0.0526723 0.0487431	0.109046 -0.0557527 0.0639262	0.0935591 -0.0526723 0.0487431	0.0953216 -0.0484171 0.0639262	0.0953216 -0.0484171 0.0639262	0.0935591 -0.0526723 0.0487431	0.0795317 -0.0484171 0.0487431	0.0953216 -0.0484171 0.0639262	0.0795317 -0.0484171 0.0487431	0.0804302 -0.0438999 0.0639262	0.0804302 -0.0438999 0.0639262	0.0795317 -0.0484171 0.0487431	0.0649438 -0.0469804 0.0487431	0.0804302 -0.0438999 0.0639262	0.0649438 -0.0469804 0.0487431	0.0649438 -0.0423746 0.0639262	0.0649438 -0.0469804 0.0487431	0.0649438 -0.0544597 0.0347502	0.051815 -0.0557528 0.0347502	0.0649438 -0.0469804 0.0487431	0.051815 -0.0557528 0.0347502	0.0503558 -0.0484172 0.0487431	0.0503558 -0.0484172 0.0487431	0.051815 -0.0557528 0.0347502	0.0391906 -0.0595824 0.0347502	0.0503558 -0.0484172 0.0487431	0.0391906 -0.0595824 0.0347502	0.0363284 -0.0526723 0.0487431	0.0363284 -0.0526723 0.0487431	0.0391906 -0.0595824 0.0347502	0.027556 -0.0658012 0.0347502	0.0363284 -0.0526723 0.0487431	0.027556 -0.0658012 0.0347502	0.0234007 -0.0595824 0.0487431	0.0234007 -0.0595824 0.0487431	0.027556 -0.0658012 0.0347502	0.0173581 -0.0741704 0.0347502	0.0234007 -0.0595824 0.0487431	0.0173581 -0.0741704 0.0347502	0.0120694 -0.0688817 0.0487431	0.0120694 -0.0688817 0.0487431	0.0173581 -0.0741704 0.0347502	0.00898895 -0.0843682 0.0347502	0.0120694 -0.0688817 0.0487431	0.00898895 -0.0843682 0.0347502	0.00277011 -0.0802129 0.0487431	0.00277011 -0.0802129 0.0487431	0.00898895 -0.0843682 0.0347502	0.0027701 -0.0960029 0.0347502	0.00277011 -0.0802129 0.0487431	0.0027701 -0.0960029 0.0347502	-0.00413991 -0.0931407 0.0487431	-0.00413991 -0.0931407 0.0487431	0.0027701 -0.0960029 0.0347502	-0.00105944 -0.108627 0.0347502	-0.00413991 -0.0931407 0.0487431	-0.00105944 -0.108627 0.0347502	-0.00839507 -0.107168 0.0487431	-0.00839507 -0.107168 0.0487431	-0.00105944 -0.108627 0.0347502	-0.00235252 -0.121756 0.0347502	-0.00839507 -0.107168 0.0487431	-0.00235252 -0.121756 0.0347502	-0.00983187 -0.121756 0.0487431	-0.00983187 -0.121756 0.0487431	-0.00235252 -0.121756 0.0347502	-0.00105944 -0.134885 0.0347502	-0.00983187 -0.121756 0.0487431	-0.00105944 -0.134885 0.0347502	-0.00839507 -0.136344 0.0487431	-0.00839507 -0.136344 0.0487431	-0.00105944 -0.134885 0.0347502	0.00277012 -0.147509 0.0347502	-0.00839507 -0.136344 0.0487431	0.00277012 -0.147509 0.0347502	-0.00413989 -0.150372 0.0487431	-0.00413989 -0.150372 0.0487431	0.00277012 -0.147509 0.0347502	0.00898898 -0.159144 0.0347502	-0.00413989 -0.150372 0.0487431	0.00898898 -0.159144 0.0347502	0.00277013 -0.163299 0.0487431	0.00277013 -0.163299 0.0487431	0.00898898 -0.159144 0.0347502	0.0173581 -0.169342 0.0347502	0.00277013 -0.163299 0.0487431	0.0173581 -0.169342 0.0347502	0.0120695 -0.17463 0.0487431	0.0120695 -0.17463 0.0487431	0.0173581 -0.169342 0.0347502	0.027556 -0.177711 0.0347502	0.0120695 -0.17463 0.0487431	0.027556 -0.177711 0.0347502	0.0234007 -0.18393 0.0487431	0.0234007 -0.18393 0.0487431	0.027556 -0.177711 0.0347502	0.0391907 -0.18393 0.0347502	0.0234007 -0.18393 0.0487431	0.0391907 -0.18393 0.0347502	0.0363285 -0.19084 0.0487431	0.0363285 -0.19084 0.0487431	0.0391907 -0.18393 0.0347502	0.051815 -0.187759 0.0347502	0.0363285 -0.19084 0.0487431	0.051815 -0.187759 0.0347502	0.0503559 -0.195095 0.0487431	0.0503559 -0.195095 0.0487431	0.051815 -0.187759 0.0347502	0.0649439 -0.189052 0.0347502	0.0503559 -0.195095 0.0487431	0.0649439 -0.189052 0.0347502	0.0649439 -0.196532 0.0487431	0.0649439 -0.196532 0.0487431	0.0649439 -0.189052 0.0347502	0.0780727 -0.187759 0.0347502	0.0649439 -0.196532 0.0487431	0.0780727 -0.187759 0.0347502	0.0795319 -0.195095 0.0487431	0.0795319 -0.195095 0.0487431	0.0780727 -0.187759 0.0347502	0.0906971 -0.18393 0.0347502	0.0795319 -0.195095 0.0487431	0.0906971 -0.18393 0.0347502	0.0935593 -0.19084 0.0487431	0.0935593 -0.19084 0.0487431	0.0906971 -0.18393 0.0347502	0.102332 -0.177711 0.0347502	0.0935593 -0.19084 0.0487431	0.102332 -0.177711 0.0347502	0.106487 -0.18393 0.0487431	0.106487 -0.18393 0.0487431	0.102332 -0.177711 0.0347502	0.11253 -0.169342 0.0347502	0.106487 -0.18393 0.0487431	0.11253 -0.169342 0.0347502	0.117818 -0.17463 0.0487431	0.117818 -0.17463 0.0487431	0.11253 -0.169342 0.0347502	0.120899 -0.159144 0.0347502	0.117818 -0.17463 0.0487431	0.120899 -0.159144 0.0347502	0.127118 -0.163299 0.0487431	0.127118 -0.163299 0.0487431	0.120899 -0.159144 0.0347502	0.127118 -0.147509 0.0347502	0.127118 -0.163299 0.0487431	0.127118 -0.147509 0.0347502	0.134028 -0.150371 0.0487431	0.134028 -0.150371 0.0487431	0.127118 -0.147509 0.0347502	0.130947 -0.134885 0.0347502	0.134028 -0.150371 0.0487431	0.130947 -0.134885 0.0347502	0.138283 -0.136344 0.0487431	0.138283 -0.136344 0.0487431	0.130947 -0.134885 0.0347502	0.13224 -0.121756 0.0347502	0.138283 -0.136344 0.0487431	0.13224 -0.121756 0.0347502	0.139719 -0.121756 0.0487431	0.139719 -0.121756 0.0487431	0.13224 -0.121756 0.0347502	0.130947 -0.108627 0.0347502	0.139719 -0.121756 0.0487431	0.130947 -0.108627 0.0347502	0.138283 -0.107168 0.0487431	0.138283 -0.107168 0.0487431	0.130947 -0.108627 0.0347502	0.127117 -0.0960028 0.0347502	0.138283 -0.107168 0.0487431	0.127117 -0.0960028 0.0347502	0.134027 -0.0931405 0.0487431	0.134027 -0.0931405 0.0487431	0.127117 -0.0960028 0.0347502	0.120899 -0.0843681 0.0347502	0.134027 -0.0931405 0.0487431	0.120899 -0.0843681 0.0347502	0.127117 -0.0802128 0.0487431	0.127117 -0.0802128 0.0487431	0.120899 -0.0843681 0.0347502	0.112529 -0.0741703 0.0347502	0.127117 -0.0802128 0.0487431	0.112529 -0.0741703 0.0347502	0.117818 -0.0688816 0.0487431	0.117818 -0.0688816 0.0487431	0.112529 -0.0741703 0.0347502	0.102332 -0.0658011 0.0347502	0.117818 -0.0688816 0.0487431	0.102332 -0.0658011 0.0347502	0.106487 -0.0595823 0.0487431	0.106487 -0.0595823 0.0487431	0.102332 -0.0658011 0.0347502	0.0906969 -0.0595823 0.0347502	0.106487 -0.0595823 0.0487431	0.0906969 -0.0595823 0.0347502	0.0935591 -0.0526723 0.0487431	0.0935591 -0.0526723 0.0487431	0.0906969 -0.0595823 0.0347502	0.0780725 -0.0557528 0.0347502	0.0935591 -0.0526723 0.0487431	0.0780725 -0.0557528 0.0347502	0.0795317 -0.0484171 0.0487431	0.0795317 -0.0484171 0.0487431	0.0780725 -0.0557528 0.0347502	0.0649438 -0.0544597 0.0347502	0.0795317 -0.0484171 0.0487431	0.0649438 -0.0544597 0.0347502	0.0649438 -0.0469804 0.0487431	0.0649438 -0.0544597 0.0347502	0.0649438 -0.0645252 0.0224853	0.0537786 -0.0656249 0.0224853	0.0649438 -0.0544597 0.0347502	0.0537786 -0.0656249 0.0224853	0.051815 -0.0557528 0.0347502	0.051815 -0.0557528 0.0347502	0.0537786 -0.0656249 0.0224853	0.0430425 -0.0688817 0.0224853	0.051815 -0.0557528 0.0347502	0.0430425 -0.0688817 0.0224853	0.0391906 -0.0595824 0.0347502	0.0391906 -0.0595824 0.0347502	0.0430425 -0.0688817 0.0224853	0.0331481 -0.0741704 0.0224853	0.0391906 -0.0595824 0.0347502	0.0331481 -0.0741704 0.0224853	0.027556 -0.0658012 0.0347502	0.027556 -0.0658012 0.0347502	0.0331481 -0.0741704 0.0224853	0.0244755 -0.0812878 0.0224853	0.027556 -0.0658012 0.0347502	0.0244755 -0.0812878 0.0224853	0.0173581 -0.0741704 0.0347502	0.0173581 -0.0741704 0.0347502	0.0244755 -0.0812878 0.0224853	0.0173581 -0.0899603 0.0224853	0.0173581 -0.0741704 0.0347502	0.0173581 -0.0899603 0.0224853	0.00898895 -0.0843682 0.0347502	0.00898895 -0.0843682 0.0347502	0.0173581 -0.0899603 0.0224853	0.0120694 -0.0998548 0.0224853	0.00898895 -0.0843682 0.0347502	0.0120694 -0.0998548 0.0224853	0.0027701 -0.0960029 0.0347502	0.0027701 -0.0960029 0.0347502	0.0120694 -0.0998548 0.0224853	0.00881266 -0.110591 0.0224853	0.0027701 -0.0960029 0.0347502	0.00881266 -0.110591 0.0224853	-0.00105944 -0.108627 0.0347502	-0.00105944 -0.108627 0.0347502	0.00881266 -0.110591 0.0224853	0.00771299 -0.121756 0.0224853	-0.00105944 -0.108627 0.0347502	0.00771299 -0.121756 0.0224853	-0.00235252 -0.121756 0.0347502	-0.00235252 -0.121756 0.0347502	0.00771299 -0.121756 0.0224853	0.00881267 -0.132921 0.0224853	-0.00235252 -0.121756 0.0347502	0.00881267 -0.132921 0.0224853	-0.00105944 -0.134885 0.0347502	-0.00105944 -0.134885 0.0347502	0.00881267 -0.132921 0.0224853	0.0120694 -0.143657 0.0224853	-0.00105944 -0.134885 0.0347502	0.0120694 -0.143657 0.0224853	0.00277012 -0.147509 0.0347502	0.00277012 -0.147509 0.0347502	0.0120694 -0.143657 0.0224853	0.0173581 -0.153552 0.0224853	0.00277012 -0.147509 0.0347502	0.0173581 -0.153552 0.0224853	0.00898898 -0.159144 0.0347502	0.00898898 -0.159144 0.0347502	0.0173581 -0.153552 0.0224853	0.0244755 -0.162224 0.0224853	0.00898898 -0.159144 0.0347502	0.0244755 -0.162224 0.0224853	0.0173581 -0.169342 0.0347502	0.0173581 -0.169342 0.0347502	0.0244755 -0.162224 0.0224853	0.0331481 -0.169342 0.0224853	0.0173581 -0.169342 0.0347502	0.0331481 -0.169342 0.0224853	0.027556 -0.177711 0.0347502	0.027556 -0.177711 0.0347502	0.0331481 -0.169342 0.0224853	0.0430426 -0.17463 0.0224853	0.027556 -0.177711 0.0347502	0.0430426 -0.17463 0.0224853	0.0391907 -0.18393 0.0347502	0.0391907 -0.18393 0.0347502	0.0430426 -0.17463 0.0224853	0.0537787 -0.177887 0.0224853	0.0391907 -0.18393 0.0347502	0.0537787 -0.177887 0.0224853	0.051815 -0.187759 0.0347502	0.051815 -0.187759 0.0347502	0.0537787 -0.177887 0.0224853	0.0649439 -0.178987 0.0224853	0.051815 -0.187759 0.0347502	0.0649439 -0.178987 0.0224853	0.0649439 -0.189052 0.0347502	0.0649439 -0.189052 0.0347502	0.0649439 -0.178987 0.0224853	0.0761091 -0.177887 0.0224853	0.0649439 -0.189052 0.0347502	0.0761091 -0.177887 0.0224853	0.0780727 -0.187759 0.0347502	0.0780727 -0.187759 0.0347502	0.0761091 -0.177887 0.0224853	0.0868452 -0.17463 0.0224853	0.0780727 -0.187759 0.0347502	0.0868452 -0.17463 0.0224853	0.0906971 -0.18393 0.0347502	0.0906971 -0.18393 0.0347502	0.0868452 -0.17463 0.0224853	0.0967396 -0.169342 0.0224853	0.0906971 -0.18393 0.0347502	0.0967396 -0.169342 0.0224853	0.102332 -0.177711 0.0347502	0.102332 -0.177711 0.0347502	0.0967396 -0.169342 0.0224853	0.105412 -0.162224 0.0224853	0.102332 -0.177711 0.0347502	0.105412 -0.162224 0.0224853	0.11253 -0.169342 0.0347502	0.11253 -0.169342 0.0347502	0.105412 -0.162224 0.0224853	0.11253 -0.153552 0.0224853	0.11253 -0.169342 0.0347502	0.11253 -0.153552 0.0224853	0.120899 -0.159144 0.0347502	0.120899 -0.159144 0.0347502	0.11253 -0.153552 0.0224853	0.117818 -0.143657 0.0224853	0.120899 -0.159144 0.0347502	0.117818 -0.143657 0.0224853	0.127118 -0.147509 0.0347502	0.127118 -0.147509 0.0347502	0.117818 -0.143657 0.0224853	0.121075 -0.132921 0.0224853	0.127118 -0.147509 0.0347502	0.121075 -0.132921 0.0224853	0.130947 -0.134885 0.0347502	0.130947 -0.134885 0.0347502	0.121075 -0.132921 0.0224853	0.122175 -0.121756 0.0224853	0.130947 -0.134885 0.0347502	0.122175 -0.121756 0.0224853	0.13224 -0.121756 0.0347502	0.13224 -0.121756 0.0347502	0.122175 -0.121756 0.0224853	0.121075 -0.110591 0.0224853	0.13224 -0.121756 0.0347502	0.121075 -0.110591 0.0224853	0.130947 -0.108627 0.0347502	0.130947 -0.108627 0.0347502	0.121075 -0.110591 0.0224853	0.117818 -0.0998547 0.0224853	0.130947 -0.108627 0.0347502	0.117818 -0.0998547 0.0224853	0.127117 -0.0960028 0.0347502	0.127117 -0.0960028 0.0347502	0.117818 -0.0998547 0.0224853	0.112529 -0.0899602 0.0224853	0.127117 -0.0960028 0.0347502	0.112529 -0.0899602 0.0224853	0.120899 -0.0843681 0.0347502	0.120899 -0.0843681 0.0347502	0.112529 -0.0899602 0.0224853	0.105412 -0.0812877 0.0224853	0.120899 -0.0843681 0.0347502	0.105412 -0.0812877 0.0224853	0.112529 -0.0741703 0.0347502	0.112529 -0.0741703 0.0347502	0.105412 -0.0812877 0.0224853	0.0967395 -0.0741703 0.0224853	0.112529 -0.0741703 0.0347502	0.0967395 -0.0741703 0.0224853	0.102332 -0.0658011 0.0347502	0.102332 -0.0658011 0.0347502	0.0967395 -0.0741703 0.0224853	0.086845 -0.0688816 0.0224853	0.102332 -0.0658011 0.0347502	0.086845 -0.0688816 0.0224853	0.0906969 -0.0595823 0.0347502	0.0906969 -0.0595823 0.0347502	0.086845 -0.0688816 0.0224853	0.0761089 -0.0656249 0.0224853	0.0906969 -0.0595823 0.0347502	0.0761089 -0.0656249 0.0224853	0.0780725 -0.0557528 0.0347502	0.0780725 -0.0557528 0.0347502	0.0761089 -0.0656249 0.0224853	0.0649438 -0.0645252 0.0224853	0.0780725 -0.0557528 0.0347502	0.0649438 -0.0645252 0.0224853	0.0649438 -0.0544597 0.0347502	0.0649438 -0.0645252 0.0224853	0.0649438 -0.0767901 0.0124198	0.0561714 -0.0776541 0.0124198	0.0649438 -0.0645252 0.0224853	0.0561714 -0.0776541 0.0124198	0.0537786 -0.0656249 0.0224853	0.0537786 -0.0656249 0.0224853	0.0561714 -0.0776541 0.0124198	0.0477361 -0.0802129 0.0124198	0.0537786 -0.0656249 0.0224853	0.0477361 -0.0802129 0.0124198	0.0430425 -0.0688817 0.0224853	0.0430425 -0.0688817 0.0224853	0.0477361 -0.0802129 0.0124198	0.0399621 -0.0843682 0.0124198	0.0430425 -0.0688817 0.0224853	0.0399621 -0.0843682 0.0124198	0.0331481 -0.0741704 0.0224853	0.0331481 -0.0741704 0.0224853	0.0399621 -0.0843682 0.0124198	0.0331481 -0.0899603 0.0124198	0.0331481 -0.0741704 0.0224853	0.0331481 -0.0899603 0.0124198	0.0244755 -0.0812878 0.0224853	0.0244755 -0.0812878 0.0224853	0.0331481 -0.0899603 0.0124198	0.027556 -0.0967743 0.0124198	0.0244755 -0.0812878 0.0224853	0.027556 -0.0967743 0.0124198	0.0173581 -0.0899603 0.0224853	0.0173581 -0.0899603 0.0224853	0.027556 -0.0967743 0.0124198	0.0234007 -0.104548 0.0124198	0.0173581 -0.0899603 0.0224853	0.0234007 -0.104548 0.0124198	0.0120694 -0.0998548 0.0224853	0.0120694 -0.0998548 0.0224853	0.0234007 -0.104548 0.0124198	0.0208419 -0.112984 0.0124198	0.0120694 -0.0998548 0.0224853	0.0208419 -0.112984 0.0124198	0.00881266 -0.110591 0.0224853	0.00881266 -0.110591 0.0224853	0.0208419 -0.112984 0.0124198	0.0199778 -0.121756 0.0124198	0.00881266 -0.110591 0.0224853	0.0199778 -0.121756 0.0124198	0.00771299 -0.121756 0.0224853	0.00771299 -0.121756 0.0224853	0.0199778 -0.121756 0.0124198	0.0208419 -0.130529 0.0124198	0.00771299 -0.121756 0.0224853	0.0208419 -0.130529 0.0124198	0.00881267 -0.132921 0.0224853	0.00881267 -0.132921 0.0224853	0.0208419 -0.130529 0.0124198	0.0234007 -0.138964 0.0124198	0.00881267 -0.132921 0.0224853	0.0234007 -0.138964 0.0124198	0.0120694 -0.143657 0.0224853	0.0120694 -0.143657 0.0224853	0.0234007 -0.138964 0.0124198	0.027556 -0.146738 0.0124198	0.0120694 -0.143657 0.0224853	0.027556 -0.146738 0.0124198	0.0173581 -0.153552 0.0224853	0.0173581 -0.153552 0.0224853	0.027556 -0.146738 0.0124198	0.0331481 -0.153552 0.0124198	0.0173581 -0.153552 0.0224853	0.0331481 -0.153552 0.0124198	0.0244755 -0.162224 0.0224853	0.0244755 -0.162224 0.0224853	0.0331481 -0.153552 0.0124198	0.0399621 -0.159144 0.0124198	0.0244755 -0.162224 0.0224853	0.0399621 -0.159144 0.0124198	0.0331481 -0.169342 0.0224853	0.0331481 -0.169342 0.0224853	0.0399621 -0.159144 0.0124198	0.0477361 -0.163299 0.0124198	0.0331481 -0.169342 0.0224853	0.0477361 -0.163299 0.0124198	0.0430426 -0.17463 0.0224853	0.0430426 -0.17463 0.0224853	0.0477361 -0.163299 0.0124198	0.0561714 -0.165858 0.0124198	0.0430426 -0.17463 0.0224853	0.0561714 -0.165858 0.0124198	0.0537787 -0.177887 0.0224853	0.0537787 -0.177887 0.0224853	0.0561714 -0.165858 0.0124198	0.0649439 -0.166722 0.0124198	0.0537787 -0.177887 0.0224853	0.0649439 -0.166722 0.0124198	0.0649439 -0.178987 0.0224853	0.0649439 -0.178987 0.0224853	0.0649439 -0.166722 0.0124198	0.0737163 -0.165858 0.0124198	0.0649439 -0.178987 0.0224853	0.0737163 -0.165858 0.0124198	0.0761091 -0.177887 0.0224853	0.0761091 -0.177887 0.0224853	0.0737163 -0.165858 0.0124198	0.0821516 -0.163299 0.0124198	0.0761091 -0.177887 0.0224853	0.0821516 -0.163299 0.0124198	0.0868452 -0.17463 0.0224853	0.0868452 -0.17463 0.0224853	0.0821516 -0.163299 0.0124198	0.0899256 -0.159144 0.0124198	0.0868452 -0.17463 0.0224853	0.0899256 -0.159144 0.0124198	0.0967396 -0.169342 0.0224853	0.0967396 -0.169342 0.0224853	0.0899256 -0.159144 0.0124198	0.0967396 -0.153552 0.0124198	0.0967396 -0.169342 0.0224853	0.0967396 -0.153552 0.0124198	0.105412 -0.162224 0.0224853	0.105412 -0.162224 0.0224853	0.0967396 -0.153552 0.0124198	0.102332 -0.146738 0.0124198	0.105412 -0.162224 0.0224853	0.102332 -0.146738 0.0124198	0.11253 -0.153552 0.0224853	0.11253 -0.153552 0.0224853	0.102332 -0.146738 0.0124198	0.106487 -0.138964 0.0124198	0.11253 -0.153552 0.0224853	0.106487 -0.138964 0.0124198	0.117818 -0.143657 0.0224853	0.117818 -0.143657 0.0224853	0.106487 -0.138964 0.0124198	0.109046 -0.130528 0.0124198	0.117818 -0.143657 0.0224853	0.109046 -0.130528 0.0124198	0.121075 -0.132921 0.0224853	0.121075 -0.132921 0.0224853	0.109046 -0.130528 0.0124198	0.10991 -0.121756 0.0124198	0.121075 -0.132921 0.0224853	0.10991 -0.121756 0.0124198	0.122175 -0.121756 0.0224853	0.122175 -0.121756 0.0224853	0.10991 -0.121756 0.0124198	0.109046 -0.112984 0.0124198	0.122175 -0.121756 0.0224853	0.109046 -0.112984 0.0124198	0.121075 -0.110591 0.0224853	0.121075 -0.110591 0.0224853	0.109046 -0.112984 0.0124198	0.106487 -0.104548 0.0124198	0.121075 -0.110591 0.0224853	0.106487 -0.104548 0.0124198	0.117818 -0.0998547 0.0224853	0.117818 -0.0998547 0.0224853	0.106487 -0.104548 0.0124198	0.102332 -0.0967742 0.0124198	0.117818 -0.0998547 0.0224853	0.102332 -0.0967742 0.0124198	0.112529 -0.0899602 0.0224853	0.112529 -0.0899602 0.0224853	0.102332 -0.0967742 0.0124198	0.0967395 -0.0899603 0.0124198	0.112529 -0.0899602 0.0224853	0.0967395 -0.0899603 0.0124198	0.105412 -0.0812877 0.0224853	0.105412 -0.0812877 0.0224853	0.0967395 -0.0899603 0.0124198	0.0899255 -0.0843682 0.0124198	0.105412 -0.0812877 0.0224853	0.0899255 -0.0843682 0.0124198	0.0967395 -0.0741703 0.0224853	0.0967395 -0.0741703 0.0224853	0.0899255 -0.0843682 0.0124198	0.0821515 -0.0802129 0.0124198	0.0967395 -0.0741703 0.0224853	0.0821515 -0.0802129 0.0124198	0.086845 -0.0688816 0.0224853	0.086845 -0.0688816 0.0224853	0.0821515 -0.0802129 0.0124198	0.0737161 -0.0776541 0.0124198	0.086845 -0.0688816 0.0224853	0.0737161 -0.0776541 0.0124198	0.0761089 -0.0656249 0.0224853	0.0761089 -0.0656249 0.0224853	0.0737161 -0.0776541 0.0124198	0.0649438 -0.0767901 0.0124198	0.0761089 -0.0656249 0.0224853	0.0649438 -0.0767901 0.0124198	0.0649438 -0.0645252 0.0224853	0.0649438 -0.0767901 0.0124198	0.0649438 -0.090783 0.00494048	0.0589013 -0.0913781 0.00494048	0.0649438 -0.0767901 0.0124198	0.0589013 -0.0913781 0.00494048	0.0561714 -0.0776541 0.0124198	0.0561714 -0.0776541 0.0124198	0.0589013 -0.0913781 0.00494048	0.0530909 -0.0931407 0.00494048	0.0561714 -0.0776541 0.0124198	0.0530909 -0.0931407 0.00494048	0.0477361 -0.0802129 0.0124198	0.0477361 -0.0802129 0.0124198	0.0530909 -0.0931407 0.00494048	0.0477361 -0.0960029 0.00494048	0.0477361 -0.0802129 0.0124198	0.0477361 -0.0960029 0.00494048	0.0399621 -0.0843682 0.0124198	0.0399621 -0.0843682 0.0124198	0.0477361 -0.0960029 0.00494048	0.0430425 -0.0998548 0.00494048	0.0399621 -0.0843682 0.0124198	0.0430425 -0.0998548 0.00494048	0.0331481 -0.0899603 0.0124198	0.0331481 -0.0899603 0.0124198	0.0430425 -0.0998548 0.00494048	0.0391906 -0.104548 0.00494048	0.0331481 -0.0899603 0.0124198	0.0391906 -0.104548 0.00494048	0.027556 -0.0967743 0.0124198	0.027556 -0.0967743 0.0124198	0.0391906 -0.104548 0.00494048	0.0363284 -0.109903 0.00494048	0.027556 -0.0967743 0.0124198	0.0363284 -0.109903 0.00494048	0.0234007 -0.104548 0.0124198	0.0234007 -0.104548 0.0124198	0.0363284 -0.109903 0.00494048	0.0345659 -0.115714 0.00494048	0.0234007 -0.104548 0.0124198	0.0345659 -0.115714 0.00494048	0.0208419 -0.112984 0.0124198	0.0208419 -0.112984 0.0124198	0.0345659 -0.115714 0.00494048	0.0339707 -0.121756 0.00494048	0.0208419 -0.112984 0.0124198	0.0339707 -0.121756 0.00494048	0.0199778 -0.121756 0.0124198	0.0199778 -0.121756 0.0124198	0.0339707 -0.121756 0.00494048	0.0345659 -0.127799 0.00494048	0.0199778 -0.121756 0.0124198	0.0345659 -0.127799 0.00494048	0.0208419 -0.130529 0.0124198	0.0208419 -0.130529 0.0124198	0.0345659 -0.127799 0.00494048	0.0363284 -0.133609 0.00494048	0.0208419 -0.130529 0.0124198	0.0363284 -0.133609 0.00494048	0.0234007 -0.138964 0.0124198	0.0234007 -0.138964 0.0124198	0.0363284 -0.133609 0.00494048	0.0391906 -0.138964 0.00494048	0.0234007 -0.138964 0.0124198	0.0391906 -0.138964 0.00494048	0.027556 -0.146738 0.0124198	0.027556 -0.146738 0.0124198	0.0391906 -0.138964 0.00494048	0.0430425 -0.143657 0.00494048	0.027556 -0.146738 0.0124198	0.0430425 -0.143657 0.00494048	0.0331481 -0.153552 0.0124198	0.0331481 -0.153552 0.0124198	0.0430425 -0.143657 0.00494048	0.0477361 -0.147509 0.00494048	0.0331481 -0.153552 0.0124198	0.0477361 -0.147509 0.00494048	0.0399621 -0.159144 0.0124198	0.0399621 -0.159144 0.0124198	0.0477361 -0.147509 0.00494048	0.0530909 -0.150371 0.00494048	0.0399621 -0.159144 0.0124198	0.0530909 -0.150371 0.00494048	0.0477361 -0.163299 0.0124198	0.0477361 -0.163299 0.0124198	0.0530909 -0.150371 0.00494048	0.0589013 -0.152134 0.00494048	0.0477361 -0.163299 0.0124198	0.0589013 -0.152134 0.00494048	0.0561714 -0.165858 0.0124198	0.0561714 -0.165858 0.0124198	0.0589013 -0.152134 0.00494048	0.0649438 -0.152729 0.00494048	0.0561714 -0.165858 0.0124198	0.0649438 -0.152729 0.00494048	0.0649439 -0.166722 0.0124198	0.0649439 -0.166722 0.0124198	0.0649438 -0.152729 0.00494048	0.0709864 -0.152134 0.00494048	0.0649439 -0.166722 0.0124198	0.0709864 -0.152134 0.00494048	0.0737163 -0.165858 0.0124198	0.0737163 -0.165858 0.0124198	0.0709864 -0.152134 0.00494048	0.0767967 -0.150371 0.00494048	0.0737163 -0.165858 0.0124198	0.0767967 -0.150371 0.00494048	0.0821516 -0.163299 0.0124198	0.0821516 -0.163299 0.0124198	0.0767967 -0.150371 0.00494048	0.0821516 -0.147509 0.00494048	0.0821516 -0.163299 0.0124198	0.0821516 -0.147509 0.00494048	0.0899256 -0.159144 0.0124198	0.0899256 -0.159144 0.0124198	0.0821516 -0.147509 0.00494048	0.0868451 -0.143657 0.00494048	0.0899256 -0.159144 0.0124198	0.0868451 -0.143657 0.00494048	0.0967396 -0.153552 0.0124198	0.0967396 -0.153552 0.0124198	0.0868451 -0.143657 0.00494048	0.090697 -0.138964 0.00494048	0.0967396 -0.153552 0.0124198	0.090697 -0.138964 0.00494048	0.102332 -0.146738 0.0124198	0.102332 -0.146738 0.0124198	0.090697 -0.138964 0.00494048	0.0935592 -0.133609 0.00494048	0.102332 -0.146738 0.0124198	0.0935592 -0.133609 0.00494048	0.106487 -0.138964 0.0124198	0.106487 -0.138964 0.0124198	0.0935592 -0.133609 0.00494048	0.0953218 -0.127799 0.00494048	0.106487 -0.138964 0.0124198	0.0953218 -0.127799 0.00494048	0.109046 -0.130528 0.0124198	0.109046 -0.130528 0.0124198	0.0953218 -0.127799 0.00494048	0.0959169 -0.121756 0.00494048	0.109046 -0.130528 0.0124198	0.0959169 -0.121756 0.00494048	0.10991 -0.121756 0.0124198	0.10991 -0.121756 0.0124198	0.0959169 -0.121756 0.00494048	0.0953218 -0.115713 0.00494048	0.10991 -0.121756 0.0124198	0.0953218 -0.115713 0.00494048	0.109046 -0.112984 0.0124198	0.109046 -0.112984 0.0124198	0.0953218 -0.115713 0.00494048	0.0935592 -0.109903 0.00494048	0.109046 -0.112984 0.0124198	0.0935592 -0.109903 0.00494048	0.106487 -0.104548 0.0124198	0.106487 -0.104548 0.0124198	0.0935592 -0.109903 0.00494048	0.090697 -0.104548 0.00494048	0.106487 -0.104548 0.0124198	0.090697 -0.104548 0.00494048	0.102332 -0.0967742 0.0124198	0.102332 -0.0967742 0.0124198	0.090697 -0.104548 0.00494048	0.086845 -0.0998547 0.00494048	0.102332 -0.0967742 0.0124198	0.086845 -0.0998547 0.00494048	0.0967395 -0.0899603 0.0124198	0.0967395 -0.0899603 0.0124198	0.086845 -0.0998547 0.00494048	0.0821515 -0.0960028 0.00494048	0.0967395 -0.0899603 0.0124198	0.0821515 -0.0960028 0.00494048	0.0899255 -0.0843682 0.0124198	0.0899255 -0.0843682 0.0124198	0.0821515 -0.0960028 0.00494048	0.0767966 -0.0931406 0.00494048	0.0899255 -0.0843682 0.0124198	0.0767966 -0.0931406 0.00494048	0.0821515 -0.0802129 0.0124198	0.0821515 -0.0802129 0.0124198	0.0767966 -0.0931406 0.00494048	0.0709863 -0.0913781 0.00494048	0.0821515 -0.0802129 0.0124198	0.0709863 -0.0913781 0.00494048	0.0737161 -0.0776541 0.0124198	0.0737161 -0.0776541 0.0124198	0.0709863 -0.0913781 0.00494048	0.0649438 -0.090783 0.00494048	0.0737161 -0.0776541 0.0124198	0.0649438 -0.090783 0.00494048	0.0649438 -0.0767901 0.0124198	0.0649438 -0.090783 0.00494048	0.0649438 -0.105966 0.000334725	0.0618634 -0.10627 0.000334725	0.0649438 -0.090783 0.00494048	0.0618634 -0.10627 0.000334725	0.0589013 -0.0913781 0.00494048	0.0589013 -0.0913781 0.00494048	0.0618634 -0.10627 0.000334725	0.0589013 -0.107168 0.000334725	0.0589013 -0.0913781 0.00494048	0.0589013 -0.107168 0.000334725	0.0530909 -0.0931407 0.00494048	0.0530909 -0.0931407 0.00494048	0.0589013 -0.107168 0.000334725	0.0561714 -0.108627 0.000334725	0.0530909 -0.0931407 0.00494048	0.0561714 -0.108627 0.000334725	0.0477361 -0.0960029 0.00494048	0.0477361 -0.0960029 0.00494048	0.0561714 -0.108627 0.000334725	0.0537786 -0.110591 0.000334725	0.0477361 -0.0960029 0.00494048	0.0537786 -0.110591 0.000334725	0.0430425 -0.0998548 0.00494048	0.0430425 -0.0998548 0.00494048	0.0537786 -0.110591 0.000334725	0.051815 -0.112984 0.000334725	0.0430425 -0.0998548 0.00494048	0.051815 -0.112984 0.000334725	0.0391906 -0.104548 0.00494048	0.0391906 -0.104548 0.00494048	0.051815 -0.112984 0.000334725	0.0503558 -0.115714 0.000334725	0.0391906 -0.104548 0.00494048	0.0503558 -0.115714 0.000334725	0.0363284 -0.109903 0.00494048	0.0363284 -0.109903 0.00494048	0.0503558 -0.115714 0.000334725	0.0494573 -0.118676 0.000334725	0.0363284 -0.109903 0.00494048	0.0494573 -0.118676 0.000334725	0.0345659 -0.115714 0.00494048	0.0345659 -0.115714 0.00494048	0.0494573 -0.118676 0.000334725	0.0491539 -0.121756 0.000334725	0.0345659 -0.115714 0.00494048	0.0491539 -0.121756 0.000334725	0.0339707 -0.121756 0.00494048	0.0339707 -0.121756 0.00494048	0.0491539 -0.121756 0.000334725	0.0494573 -0.124837 0.000334725	0.0339707 -0.121756 0.00494048	0.0494573 -0.124837 0.000334725	0.0345659 -0.127799 0.00494048	0.0345659 -0.127799 0.00494048	0.0494573 -0.124837 0.000334725	0.0503558 -0.127799 0.000334725	0.0345659 -0.127799 0.00494048	0.0503558 -0.127799 0.000334725	0.0363284 -0.133609 0.00494048	0.0363284 -0.133609 0.00494048	0.0503558 -0.127799 0.000334725	0.051815 -0.130528 0.000334725	0.0363284 -0.133609 0.00494048	0.051815 -0.130528 0.000334725	0.0391906 -0.138964 0.00494048	0.0391906 -0.138964 0.00494048	0.051815 -0.130528 0.000334725	0.0537787 -0.132921 0.000334725	0.0391906 -0.138964 0.00494048	0.0537787 -0.132921 0.000334725	0.0430425 -0.143657 0.00494048	0.0430425 -0.143657 0.00494048	0.0537787 -0.132921 0.000334725	0.0561714 -0.134885 0.000334725	0.0430425 -0.143657 0.00494048	0.0561714 -0.134885 0.000334725	0.0477361 -0.147509 0.00494048	0.0477361 -0.147509 0.00494048	0.0561714 -0.134885 0.000334725	0.0589013 -0.136344 0.000334725	0.0477361 -0.147509 0.00494048	0.0589013 -0.136344 0.000334725	0.0530909 -0.150371 0.00494048	0.0530909 -0.150371 0.00494048	0.0589013 -0.136344 0.000334725	0.0618634 -0.137243 0.000334725	0.0530909 -0.150371 0.00494048	0.0618634 -0.137243 0.000334725	0.0589013 -0.152134 0.00494048	0.0589013 -0.152134 0.00494048	0.0618634 -0.137243 0.000334725	0.0649438 -0.137546 0.000334725	0.0589013 -0.152134 0.00494048	0.0649438 -0.137546 0.000334725	0.0649438 -0.152729 0.00494048	0.0649438 -0.152729 0.00494048	0.0649438 -0.137546 0.000334725	0.0680243 -0.137243 0.000334725	0.0649438 -0.152729 0.00494048	0.0680243 -0.137243 0.000334725	0.0709864 -0.152134 0.00494048	0.0709864 -0.152134 0.00494048	0.0680243 -0.137243 0.000334725	0.0709864 -0.136344 0.000334725	0.0709864 -0.152134 0.00494048	0.0709864 -0.136344 0.000334725	0.0767967 -0.150371 0.00494048	0.0767967 -0.150371 0.00494048	0.0709864 -0.136344 0.000334725	0.0737162 -0.134885 0.000334725	0.0767967 -0.150371 0.00494048	0.0737162 -0.134885 0.000334725	0.0821516 -0.147509 0.00494048	0.0821516 -0.147509 0.00494048	0.0737162 -0.134885 0.000334725	0.076109 -0.132921 0.000334725	0.0821516 -0.147509 0.00494048	0.076109 -0.132921 0.000334725	0.0868451 -0.143657 0.00494048	0.0868451 -0.143657 0.00494048	0.076109 -0.132921 0.000334725	0.0780727 -0.130528 0.000334725	0.0868451 -0.143657 0.00494048	0.0780727 -0.130528 0.000334725	0.090697 -0.138964 0.00494048	0.090697 -0.138964 0.00494048	0.0780727 -0.130528 0.000334725	0.0795318 -0.127799 0.000334725	0.090697 -0.138964 0.00494048	0.0795318 -0.127799 0.000334725	0.0935592 -0.133609 0.00494048	0.0935592 -0.133609 0.00494048	0.0795318 -0.127799 0.000334725	0.0804303 -0.124837 0.000334725	0.0935592 -0.133609 0.00494048	0.0804303 -0.124837 0.000334725	0.0953218 -0.127799 0.00494048	0.0953218 -0.127799 0.00494048	0.0804303 -0.124837 0.000334725	0.0807337 -0.121756 0.000334725	0.0953218 -0.127799 0.00494048	0.0807337 -0.121756 0.000334725	0.0959169 -0.121756 0.00494048	0.0959169 -0.121756 0.00494048	0.0807337 -0.121756 0.000334725	0.0804303 -0.118676 0.000334725	0.0959169 -0.121756 0.00494048	0.0804303 -0.118676 0.000334725	0.0953218 -0.115713 0.00494048	0.0953218 -0.115713 0.00494048	0.0804303 -0.118676 0.000334725	0.0795318 -0.115713 0.000334725	0.0953218 -0.115713 0.00494048	0.0795318 -0.115713 0.000334725	0.0935592 -0.109903 0.00494048	0.0935592 -0.109903 0.00494048	0.0795318 -0.115713 0.000334725	0.0780726 -0.112984 0.000334725	0.0935592 -0.109903 0.00494048	0.0780726 -0.112984 0.000334725	0.090697 -0.104548 0.00494048	0.090697 -0.104548 0.00494048	0.0780726 -0.112984 0.000334725	0.076109 -0.110591 0.000334725	0.090697 -0.104548 0.00494048	0.076109 -0.110591 0.000334725	0.086845 -0.0998547 0.00494048	0.086845 -0.0998547 0.00494048	0.076109 -0.110591 0.000334725	0.0737162 -0.108627 0.000334725	0.086845 -0.0998547 0.00494048	0.0737162 -0.108627 0.000334725	0.0821515 -0.0960028 0.00494048	0.0821515 -0.0960028 0.00494048	0.0737162 -0.108627 0.000334725	0.0709863 -0.107168 0.000334725	0.0821515 -0.0960028 0.00494048	0.0709863 -0.107168 0.000334725	0.0767966 -0.0931406 0.00494048	0.0767966 -0.0931406 0.00494048	0.0709863 -0.107168 0.000334725	0.0680242 -0.10627 0.000334725	0.0767966 -0.0931406 0.00494048	0.0680242 -0.10627 0.000334725	0.0709863 -0.0913781 0.00494048	0.0709863 -0.0913781 0.00494048	0.0680242 -0.10627 0.000334725	0.0649438 -0.105966 0.000334725	0.0709863 -0.0913781 0.00494048	0.0649438 -0.105966 0.000334725	0.0649438 -0.090783 0.00494048	0.0649438 -0.121756 -0.00122045	0.0618634 -0.10627 0.000334725	0.0649438 -0.105966 0.000334725	0.0649438 -0.121756 -0.00122045	0.0589013 -0.107168 0.000334725	0.0618634 -0.10627 0.000334725	0.0649438 -0.121756 -0.00122045	0.0561714 -0.108627 0.000334725	0.0589013 -0.107168 0.000334725	0.0649438 -0.121756 -0.00122045	0.0537786 -0.110591 0.000334725	0.0561714 -0.108627 0.000334725	0.0649438 -0.121756 -0.00122045	0.051815 -0.112984 0.000334725	0.0537786 -0.110591 0.000334725	0.0649438 -0.121756 -0.00122045	0.0503558 -0.115714 0.000334725	0.051815 -0.112984 0.000334725	0.0649438 -0.121756 -0.00122045	0.0494573 -0.118676 0.000334725	0.0503558 -0.115714 0.000334725	0.0649438 -0.121756 -0.00122045	0.0491539 -0.121756 0.000334725	0.0494573 -0.118676 0.000334725	0.0649438 -0.121756 -0.00122045	0.0494573 -0.124837 0.000334725	0.0491539 -0.121756 0.000334725	0.0649438 -0.121756 -0.00122045	0.0503558 -0.127799 0.000334725	0.0494573 -0.124837 0.000334725	0.0649438 -0.121756 -0.00122045	0.051815 -0.130528 0.000334725	0.0503558 -0.127799 0.000334725	0.0649438 -0.121756 -0.00122045	0.0537787 -0.132921 0.000334725	0.051815 -0.130528 0.000334725	0.0649438 -0.121756 -0.00122045	0.0561714 -0.134885 0.000334725	0.0537787 -0.132921 0.000334725	0.0649438 -0.121756 -0.00122045	0.0589013 -0.136344 0.000334725	0.0561714 -0.134885 0.000334725	0.0649438 -0.121756 -0.00122045	0.0618634 -0.137243 0.000334725	0.0589013 -0.136344 0.000334725	0.0649438 -0.121756 -0.00122045	0.0649438 -0.137546 0.000334725	0.0618634 -0.137243 0.000334725	0.0649438 -0.121756 -0.00122045	0.0680243 -0.137243 0.000334725	0.0649438 -0.137546 0.000334725	0.0649438 -0.121756 -0.00122045	0.0709864 -0.136344 0.000334725	0.0680243 -0.137243 0.000334725	0.0649438 -0.121756 -0.00122045	0.0737162 -0.134885 0.000334725	0.0709864 -0.136344 0.000334725	0.0649438 -0.121756 -0.00122045	0.076109 -0.132921 0.000334725	0.0737162 -0.134885 0.000334725	0.0649438 -0.121756 -0.00122045	0.0780727 -0.130528 0.000334725	0.076109 -0.132921 0.000334725	0.0649438 -0.121756 -0.00122045	0.0795318 -0.127799 0.000334725	0.0780727 -0.130528 0.000334725	0.0649438 -0.121756 -0.00122045	0.0804303 -0.124837 0.000334725	0.0795318 -0.127799 0.000334725	0.0649438 -0.121756 -0.00122045	0.0807337 -0.121756 0.000334725	0.0804303 -0.124837 0.000334725	0.0649438 -0.121756 -0.00122045	0.0804303 -0.118676 0.000334725	0.0807337 -0.121756 0.000334725	0.0649438 -0.121756 -0.00122045	0.0795318 -0.115713 0.000334725	0.0804303 -0.118676 0.000334725	0.0649438 -0.121756 -0.00122045	0.0780726 -0.112984 0.000334725	0.0795318 -0.115713 0.000334725	0.0649438 -0.121756 -0.00122045	0.076109 -0.110591 0.000334725	0.0780726 -0.112984 0.000334725	0.0649438 -0.121756 -0.00122045	0.0737162 -0.108627 0.000334725	0.076109 -0.110591 0.000334725	0.0649438 -0.121756 -0.00122045	0.0709863 -0.107168 0.000334725	0.0737162 -0.108627 0.000334725	0.0649438 -0.121756 -0.00122045	0.0680242 -0.10627 0.000334725	0.0709863 -0.107168 0.000334725	0.0649438 -0.121756 -0.00122045	0.0649438 -0.105966 0.000334725	0.0680242 -0.10627 0.000334725]\n" to:out_file
		-- opacity maps
		if classof mat == PBRT_Matte or classof mat == PBRT_Plastic or classof mat == PBRT_Translucent then (
			
				if mat.alphaAmount < 1.0 and mat.alphaEnable == false then (
						format "\t \"float alpha\" [%]\n" mat.alphaAmount  to:out_file	
					)
				if mat.alphaEnable and mat.alphaImage != undefined then (
						format "\t \"texture alpha\" \"%\"\n" mat.alphaImage.name to:out_file	
					)									
					
			)
		format "\"float uv\" [0.0 1.0	0.0 0.9375	0.03125 0.9375	0.03125 1.0	0.03125 0.9375	0.0625 0.9375	0.0625 1.0	0.0625 0.9375	0.09375 0.9375	0.09375 1.0	0.09375 0.9375	0.125 0.9375	0.125 1.0	0.125 0.9375	0.15625 0.9375	0.15625 1.0	0.15625 0.9375	0.1875 0.9375	0.1875 1.0	0.1875 0.9375	0.21875 0.9375	0.21875 1.0	0.21875 0.9375	0.25 0.9375	0.25 1.0	0.25 0.9375	0.28125 0.9375	0.28125 1.0	0.28125 0.9375	0.3125 0.9375	0.3125 1.0	0.3125 0.9375	0.34375 0.9375	0.34375 1.0	0.34375 0.9375	0.375 0.9375	0.375 1.0	0.375 0.9375	0.40625 0.9375	0.40625 1.0	0.40625 0.9375	0.4375 0.9375	0.4375 1.0	0.4375 0.9375	0.46875 0.9375	0.46875 1.0	0.46875 0.9375	0.5 0.9375	0.5 1.0	0.5 0.9375	0.53125 0.9375	0.53125 1.0	0.53125 0.9375	0.5625 0.9375	0.5625 1.0	0.5625 0.9375	0.59375 0.9375	0.59375 1.0	0.59375 0.9375	0.625 0.9375	0.625 1.0	0.625 0.9375	0.65625 0.9375	0.65625 1.0	0.65625 0.9375	0.6875 0.9375	0.6875 1.0	0.6875 0.9375	0.71875 0.9375	0.71875 1.0	0.71875 0.9375	0.75 0.9375	0.75 1.0	0.75 0.9375	0.78125 0.9375	0.78125 1.0	0.78125 0.9375	0.8125 0.9375	0.8125 1.0	0.8125 0.9375	0.84375 0.9375	0.84375 1.0	0.84375 0.9375	0.875 0.9375	0.875 1.0	0.875 0.9375	0.90625 0.9375	0.90625 1.0	0.90625 0.9375	0.9375 0.9375	0.9375 1.0	0.9375 0.9375	0.96875 0.9375	0.96875 1.0	0.96875 0.9375	1.0 0.9375	0.0 0.9375	0.0 0.875	0.03125 0.875	0.0 0.9375	0.03125 0.875	0.03125 0.9375	0.03125 0.9375	0.03125 0.875	0.0625 0.875	0.03125 0.9375	0.0625 0.875	0.0625 0.9375	0.0625 0.9375	0.0625 0.875	0.09375 0.875	0.0625 0.9375	0.09375 0.875	0.09375 0.9375	0.09375 0.9375	0.09375 0.875	0.125 0.875	0.09375 0.9375	0.125 0.875	0.125 0.9375	0.125 0.9375	0.125 0.875	0.15625 0.875	0.125 0.9375	0.15625 0.875	0.15625 0.9375	0.15625 0.9375	0.15625 0.875	0.1875 0.875	0.15625 0.9375	0.1875 0.875	0.1875 0.9375	0.1875 0.9375	0.1875 0.875	0.21875 0.875	0.1875 0.9375	0.21875 0.875	0.21875 0.9375	0.21875 0.9375	0.21875 0.875	0.25 0.875	0.21875 0.9375	0.25 0.875	0.25 0.9375	0.25 0.9375	0.25 0.875	0.28125 0.875	0.25 0.9375	0.28125 0.875	0.28125 0.9375	0.28125 0.9375	0.28125 0.875	0.3125 0.875	0.28125 0.9375	0.3125 0.875	0.3125 0.9375	0.3125 0.9375	0.3125 0.875	0.34375 0.875	0.3125 0.9375	0.34375 0.875	0.34375 0.9375	0.34375 0.9375	0.34375 0.875	0.375 0.875	0.34375 0.9375	0.375 0.875	0.375 0.9375	0.375 0.9375	0.375 0.875	0.40625 0.875	0.375 0.9375	0.40625 0.875	0.40625 0.9375	0.40625 0.9375	0.40625 0.875	0.4375 0.875	0.40625 0.9375	0.4375 0.875	0.4375 0.9375	0.4375 0.9375	0.4375 0.875	0.46875 0.875	0.4375 0.9375	0.46875 0.875	0.46875 0.9375	0.46875 0.9375	0.46875 0.875	0.5 0.875	0.46875 0.9375	0.5 0.875	0.5 0.9375	0.5 0.9375	0.5 0.875	0.53125 0.875	0.5 0.9375	0.53125 0.875	0.53125 0.9375	0.53125 0.9375	0.53125 0.875	0.5625 0.875	0.53125 0.9375	0.5625 0.875	0.5625 0.9375	0.5625 0.9375	0.5625 0.875	0.59375 0.875	0.5625 0.9375	0.59375 0.875	0.59375 0.9375	0.59375 0.9375	0.59375 0.875	0.625 0.875	0.59375 0.9375	0.625 0.875	0.625 0.9375	0.625 0.9375	0.625 0.875	0.65625 0.875	0.625 0.9375	0.65625 0.875	0.65625 0.9375	0.65625 0.9375	0.65625 0.875	0.6875 0.875	0.65625 0.9375	0.6875 0.875	0.6875 0.9375	0.6875 0.9375	0.6875 0.875	0.71875 0.875	0.6875 0.9375	0.71875 0.875	0.71875 0.9375	0.71875 0.9375	0.71875 0.875	0.75 0.875	0.71875 0.9375	0.75 0.875	0.75 0.9375	0.75 0.9375	0.75 0.875	0.78125 0.875	0.75 0.9375	0.78125 0.875	0.78125 0.9375	0.78125 0.9375	0.78125 0.875	0.8125 0.875	0.78125 0.9375	0.8125 0.875	0.8125 0.9375	0.8125 0.9375	0.8125 0.875	0.84375 0.875	0.8125 0.9375	0.84375 0.875	0.84375 0.9375	0.84375 0.9375	0.84375 0.875	0.875 0.875	0.84375 0.9375	0.875 0.875	0.875 0.9375	0.875 0.9375	0.875 0.875	0.90625 0.875	0.875 0.9375	0.90625 0.875	0.90625 0.9375	0.90625 0.9375	0.90625 0.875	0.9375 0.875	0.90625 0.9375	0.9375 0.875	0.9375 0.9375	0.9375 0.9375	0.9375 0.875	0.96875 0.875	0.9375 0.9375	0.96875 0.875	0.96875 0.9375	0.96875 0.9375	0.96875 0.875	1.0 0.875	0.96875 0.9375	1.0 0.875	1.0 0.9375	0.0 0.875	0.0 0.8125	0.03125 0.8125	0.0 0.875	0.03125 0.8125	0.03125 0.875	0.03125 0.875	0.03125 0.8125	0.0625 0.8125	0.03125 0.875	0.0625 0.8125	0.0625 0.875	0.0625 0.875	0.0625 0.8125	0.09375 0.8125	0.0625 0.875	0.09375 0.8125	0.09375 0.875	0.09375 0.875	0.09375 0.8125	0.125 0.8125	0.09375 0.875	0.125 0.8125	0.125 0.875	0.125 0.875	0.125 0.8125	0.15625 0.8125	0.125 0.875	0.15625 0.8125	0.15625 0.875	0.15625 0.875	0.15625 0.8125	0.1875 0.8125	0.15625 0.875	0.1875 0.8125	0.1875 0.875	0.1875 0.875	0.1875 0.8125	0.21875 0.8125	0.1875 0.875	0.21875 0.8125	0.21875 0.875	0.21875 0.875	0.21875 0.8125	0.25 0.8125	0.21875 0.875	0.25 0.8125	0.25 0.875	0.25 0.875	0.25 0.8125	0.28125 0.8125	0.25 0.875	0.28125 0.8125	0.28125 0.875	0.28125 0.875	0.28125 0.8125	0.3125 0.8125	0.28125 0.875	0.3125 0.8125	0.3125 0.875	0.3125 0.875	0.3125 0.8125	0.34375 0.8125	0.3125 0.875	0.34375 0.8125	0.34375 0.875	0.34375 0.875	0.34375 0.8125	0.375 0.8125	0.34375 0.875	0.375 0.8125	0.375 0.875	0.375 0.875	0.375 0.8125	0.40625 0.8125	0.375 0.875	0.40625 0.8125	0.40625 0.875	0.40625 0.875	0.40625 0.8125	0.4375 0.8125	0.40625 0.875	0.4375 0.8125	0.4375 0.875	0.4375 0.875	0.4375 0.8125	0.46875 0.8125	0.4375 0.875	0.46875 0.8125	0.46875 0.875	0.46875 0.875	0.46875 0.8125	0.5 0.8125	0.46875 0.875	0.5 0.8125	0.5 0.875	0.5 0.875	0.5 0.8125	0.53125 0.8125	0.5 0.875	0.53125 0.8125	0.53125 0.875	0.53125 0.875	0.53125 0.8125	0.5625 0.8125	0.53125 0.875	0.5625 0.8125	0.5625 0.875	0.5625 0.875	0.5625 0.8125	0.59375 0.8125	0.5625 0.875	0.59375 0.8125	0.59375 0.875	0.59375 0.875	0.59375 0.8125	0.625 0.8125	0.59375 0.875	0.625 0.8125	0.625 0.875	0.625 0.875	0.625 0.8125	0.65625 0.8125	0.625 0.875	0.65625 0.8125	0.65625 0.875	0.65625 0.875	0.65625 0.8125	0.6875 0.8125	0.65625 0.875	0.6875 0.8125	0.6875 0.875	0.6875 0.875	0.6875 0.8125	0.71875 0.8125	0.6875 0.875	0.71875 0.8125	0.71875 0.875	0.71875 0.875	0.71875 0.8125	0.75 0.8125	0.71875 0.875	0.75 0.8125	0.75 0.875	0.75 0.875	0.75 0.8125	0.78125 0.8125	0.75 0.875	0.78125 0.8125	0.78125 0.875	0.78125 0.875	0.78125 0.8125	0.8125 0.8125	0.78125 0.875	0.8125 0.8125	0.8125 0.875	0.8125 0.875	0.8125 0.8125	0.84375 0.8125	0.8125 0.875	0.84375 0.8125	0.84375 0.875	0.84375 0.875	0.84375 0.8125	0.875 0.8125	0.84375 0.875	0.875 0.8125	0.875 0.875	0.875 0.875	0.875 0.8125	0.90625 0.8125	0.875 0.875	0.90625 0.8125	0.90625 0.875	0.90625 0.875	0.90625 0.8125	0.9375 0.8125	0.90625 0.875	0.9375 0.8125	0.9375 0.875	0.9375 0.875	0.9375 0.8125	0.96875 0.8125	0.9375 0.875	0.96875 0.8125	0.96875 0.875	0.96875 0.875	0.96875 0.8125	1.0 0.8125	0.96875 0.875	1.0 0.8125	1.0 0.875	0.0 0.8125	0.0 0.75	0.03125 0.75	0.0 0.8125	0.03125 0.75	0.03125 0.8125	0.03125 0.8125	0.03125 0.75	0.0625 0.75	0.03125 0.8125	0.0625 0.75	0.0625 0.8125	0.0625 0.8125	0.0625 0.75	0.09375 0.75	0.0625 0.8125	0.09375 0.75	0.09375 0.8125	0.09375 0.8125	0.09375 0.75	0.125 0.75	0.09375 0.8125	0.125 0.75	0.125 0.8125	0.125 0.8125	0.125 0.75	0.15625 0.75	0.125 0.8125	0.15625 0.75	0.15625 0.8125	0.15625 0.8125	0.15625 0.75	0.1875 0.75	0.15625 0.8125	0.1875 0.75	0.1875 0.8125	0.1875 0.8125	0.1875 0.75	0.21875 0.75	0.1875 0.8125	0.21875 0.75	0.21875 0.8125	0.21875 0.8125	0.21875 0.75	0.25 0.75	0.21875 0.8125	0.25 0.75	0.25 0.8125	0.25 0.8125	0.25 0.75	0.28125 0.75	0.25 0.8125	0.28125 0.75	0.28125 0.8125	0.28125 0.8125	0.28125 0.75	0.3125 0.75	0.28125 0.8125	0.3125 0.75	0.3125 0.8125	0.3125 0.8125	0.3125 0.75	0.34375 0.75	0.3125 0.8125	0.34375 0.75	0.34375 0.8125	0.34375 0.8125	0.34375 0.75	0.375 0.75	0.34375 0.8125	0.375 0.75	0.375 0.8125	0.375 0.8125	0.375 0.75	0.40625 0.75	0.375 0.8125	0.40625 0.75	0.40625 0.8125	0.40625 0.8125	0.40625 0.75	0.4375 0.75	0.40625 0.8125	0.4375 0.75	0.4375 0.8125	0.4375 0.8125	0.4375 0.75	0.46875 0.75	0.4375 0.8125	0.46875 0.75	0.46875 0.8125	0.46875 0.8125	0.46875 0.75	0.5 0.75	0.46875 0.8125	0.5 0.75	0.5 0.8125	0.5 0.8125	0.5 0.75	0.53125 0.75	0.5 0.8125	0.53125 0.75	0.53125 0.8125	0.53125 0.8125	0.53125 0.75	0.5625 0.75	0.53125 0.8125	0.5625 0.75	0.5625 0.8125	0.5625 0.8125	0.5625 0.75	0.59375 0.75	0.5625 0.8125	0.59375 0.75	0.59375 0.8125	0.59375 0.8125	0.59375 0.75	0.625 0.75	0.59375 0.8125	0.625 0.75	0.625 0.8125	0.625 0.8125	0.625 0.75	0.65625 0.75	0.625 0.8125	0.65625 0.75	0.65625 0.8125	0.65625 0.8125	0.65625 0.75	0.6875 0.75	0.65625 0.8125	0.6875 0.75	0.6875 0.8125	0.6875 0.8125	0.6875 0.75	0.71875 0.75	0.6875 0.8125	0.71875 0.75	0.71875 0.8125	0.71875 0.8125	0.71875 0.75	0.75 0.75	0.71875 0.8125	0.75 0.75	0.75 0.8125	0.75 0.8125	0.75 0.75	0.78125 0.75	0.75 0.8125	0.78125 0.75	0.78125 0.8125	0.78125 0.8125	0.78125 0.75	0.8125 0.75	0.78125 0.8125	0.8125 0.75	0.8125 0.8125	0.8125 0.8125	0.8125 0.75	0.84375 0.75	0.8125 0.8125	0.84375 0.75	0.84375 0.8125	0.84375 0.8125	0.84375 0.75	0.875 0.75	0.84375 0.8125	0.875 0.75	0.875 0.8125	0.875 0.8125	0.875 0.75	0.90625 0.75	0.875 0.8125	0.90625 0.75	0.90625 0.8125	0.90625 0.8125	0.90625 0.75	0.9375 0.75	0.90625 0.8125	0.9375 0.75	0.9375 0.8125	0.9375 0.8125	0.9375 0.75	0.96875 0.75	0.9375 0.8125	0.96875 0.75	0.96875 0.8125	0.96875 0.8125	0.96875 0.75	1.0 0.75	0.96875 0.8125	1.0 0.75	1.0 0.8125	0.0 0.75	0.0 0.6875	0.03125 0.6875	0.0 0.75	0.03125 0.6875	0.03125 0.75	0.03125 0.75	0.03125 0.6875	0.0625 0.6875	0.03125 0.75	0.0625 0.6875	0.0625 0.75	0.0625 0.75	0.0625 0.6875	0.09375 0.6875	0.0625 0.75	0.09375 0.6875	0.09375 0.75	0.09375 0.75	0.09375 0.6875	0.125 0.6875	0.09375 0.75	0.125 0.6875	0.125 0.75	0.125 0.75	0.125 0.6875	0.15625 0.6875	0.125 0.75	0.15625 0.6875	0.15625 0.75	0.15625 0.75	0.15625 0.6875	0.1875 0.6875	0.15625 0.75	0.1875 0.6875	0.1875 0.75	0.1875 0.75	0.1875 0.6875	0.21875 0.6875	0.1875 0.75	0.21875 0.6875	0.21875 0.75	0.21875 0.75	0.21875 0.6875	0.25 0.6875	0.21875 0.75	0.25 0.6875	0.25 0.75	0.25 0.75	0.25 0.6875	0.28125 0.6875	0.25 0.75	0.28125 0.6875	0.28125 0.75	0.28125 0.75	0.28125 0.6875	0.3125 0.6875	0.28125 0.75	0.3125 0.6875	0.3125 0.75	0.3125 0.75	0.3125 0.6875	0.34375 0.6875	0.3125 0.75	0.34375 0.6875	0.34375 0.75	0.34375 0.75	0.34375 0.6875	0.375 0.6875	0.34375 0.75	0.375 0.6875	0.375 0.75	0.375 0.75	0.375 0.6875	0.40625 0.6875	0.375 0.75	0.40625 0.6875	0.40625 0.75	0.40625 0.75	0.40625 0.6875	0.4375 0.6875	0.40625 0.75	0.4375 0.6875	0.4375 0.75	0.4375 0.75	0.4375 0.6875	0.46875 0.6875	0.4375 0.75	0.46875 0.6875	0.46875 0.75	0.46875 0.75	0.46875 0.6875	0.5 0.6875	0.46875 0.75	0.5 0.6875	0.5 0.75	0.5 0.75	0.5 0.6875	0.53125 0.6875	0.5 0.75	0.53125 0.6875	0.53125 0.75	0.53125 0.75	0.53125 0.6875	0.5625 0.6875	0.53125 0.75	0.5625 0.6875	0.5625 0.75	0.5625 0.75	0.5625 0.6875	0.59375 0.6875	0.5625 0.75	0.59375 0.6875	0.59375 0.75	0.59375 0.75	0.59375 0.6875	0.625 0.6875	0.59375 0.75	0.625 0.6875	0.625 0.75	0.625 0.75	0.625 0.6875	0.65625 0.6875	0.625 0.75	0.65625 0.6875	0.65625 0.75	0.65625 0.75	0.65625 0.6875	0.6875 0.6875	0.65625 0.75	0.6875 0.6875	0.6875 0.75	0.6875 0.75	0.6875 0.6875	0.71875 0.6875	0.6875 0.75	0.71875 0.6875	0.71875 0.75	0.71875 0.75	0.71875 0.6875	0.75 0.6875	0.71875 0.75	0.75 0.6875	0.75 0.75	0.75 0.75	0.75 0.6875	0.78125 0.6875	0.75 0.75	0.78125 0.6875	0.78125 0.75	0.78125 0.75	0.78125 0.6875	0.8125 0.6875	0.78125 0.75	0.8125 0.6875	0.8125 0.75	0.8125 0.75	0.8125 0.6875	0.84375 0.6875	0.8125 0.75	0.84375 0.6875	0.84375 0.75	0.84375 0.75	0.84375 0.6875	0.875 0.6875	0.84375 0.75	0.875 0.6875	0.875 0.75	0.875 0.75	0.875 0.6875	0.90625 0.6875	0.875 0.75	0.90625 0.6875	0.90625 0.75	0.90625 0.75	0.90625 0.6875	0.9375 0.6875	0.90625 0.75	0.9375 0.6875	0.9375 0.75	0.9375 0.75	0.9375 0.6875	0.96875 0.6875	0.9375 0.75	0.96875 0.6875	0.96875 0.75	0.96875 0.75	0.96875 0.6875	1.0 0.6875	0.96875 0.75	1.0 0.6875	1.0 0.75	0.0 0.6875	0.0 0.625	0.03125 0.625	0.0 0.6875	0.03125 0.625	0.03125 0.6875	0.03125 0.6875	0.03125 0.625	0.0625 0.625	0.03125 0.6875	0.0625 0.625	0.0625 0.6875	0.0625 0.6875	0.0625 0.625	0.09375 0.625	0.0625 0.6875	0.09375 0.625	0.09375 0.6875	0.09375 0.6875	0.09375 0.625	0.125 0.625	0.09375 0.6875	0.125 0.625	0.125 0.6875	0.125 0.6875	0.125 0.625	0.15625 0.625	0.125 0.6875	0.15625 0.625	0.15625 0.6875	0.15625 0.6875	0.15625 0.625	0.1875 0.625	0.15625 0.6875	0.1875 0.625	0.1875 0.6875	0.1875 0.6875	0.1875 0.625	0.21875 0.625	0.1875 0.6875	0.21875 0.625	0.21875 0.6875	0.21875 0.6875	0.21875 0.625	0.25 0.625	0.21875 0.6875	0.25 0.625	0.25 0.6875	0.25 0.6875	0.25 0.625	0.28125 0.625	0.25 0.6875	0.28125 0.625	0.28125 0.6875	0.28125 0.6875	0.28125 0.625	0.3125 0.625	0.28125 0.6875	0.3125 0.625	0.3125 0.6875	0.3125 0.6875	0.3125 0.625	0.34375 0.625	0.3125 0.6875	0.34375 0.625	0.34375 0.6875	0.34375 0.6875	0.34375 0.625	0.375 0.625	0.34375 0.6875	0.375 0.625	0.375 0.6875	0.375 0.6875	0.375 0.625	0.40625 0.625	0.375 0.6875	0.40625 0.625	0.40625 0.6875	0.40625 0.6875	0.40625 0.625	0.4375 0.625	0.40625 0.6875	0.4375 0.625	0.4375 0.6875	0.4375 0.6875	0.4375 0.625	0.46875 0.625	0.4375 0.6875	0.46875 0.625	0.46875 0.6875	0.46875 0.6875	0.46875 0.625	0.5 0.625	0.46875 0.6875	0.5 0.625	0.5 0.6875	0.5 0.6875	0.5 0.625	0.53125 0.625	0.5 0.6875	0.53125 0.625	0.53125 0.6875	0.53125 0.6875	0.53125 0.625	0.5625 0.625	0.53125 0.6875	0.5625 0.625	0.5625 0.6875	0.5625 0.6875	0.5625 0.625	0.59375 0.625	0.5625 0.6875	0.59375 0.625	0.59375 0.6875	0.59375 0.6875	0.59375 0.625	0.625 0.625	0.59375 0.6875	0.625 0.625	0.625 0.6875	0.625 0.6875	0.625 0.625	0.65625 0.625	0.625 0.6875	0.65625 0.625	0.65625 0.6875	0.65625 0.6875	0.65625 0.625	0.6875 0.625	0.65625 0.6875	0.6875 0.625	0.6875 0.6875	0.6875 0.6875	0.6875 0.625	0.71875 0.625	0.6875 0.6875	0.71875 0.625	0.71875 0.6875	0.71875 0.6875	0.71875 0.625	0.75 0.625	0.71875 0.6875	0.75 0.625	0.75 0.6875	0.75 0.6875	0.75 0.625	0.78125 0.625	0.75 0.6875	0.78125 0.625	0.78125 0.6875	0.78125 0.6875	0.78125 0.625	0.8125 0.625	0.78125 0.6875	0.8125 0.625	0.8125 0.6875	0.8125 0.6875	0.8125 0.625	0.84375 0.625	0.8125 0.6875	0.84375 0.625	0.84375 0.6875	0.84375 0.6875	0.84375 0.625	0.875 0.625	0.84375 0.6875	0.875 0.625	0.875 0.6875	0.875 0.6875	0.875 0.625	0.90625 0.625	0.875 0.6875	0.90625 0.625	0.90625 0.6875	0.90625 0.6875	0.90625 0.625	0.9375 0.625	0.90625 0.6875	0.9375 0.625	0.9375 0.6875	0.9375 0.6875	0.9375 0.625	0.96875 0.625	0.9375 0.6875	0.96875 0.625	0.96875 0.6875	0.96875 0.6875	0.96875 0.625	1.0 0.625	0.96875 0.6875	1.0 0.625	1.0 0.6875	0.0 0.625	0.0 0.5625	0.03125 0.5625	0.0 0.625	0.03125 0.5625	0.03125 0.625	0.03125 0.625	0.03125 0.5625	0.0625 0.5625	0.03125 0.625	0.0625 0.5625	0.0625 0.625	0.0625 0.625	0.0625 0.5625	0.09375 0.5625	0.0625 0.625	0.09375 0.5625	0.09375 0.625	0.09375 0.625	0.09375 0.5625	0.125 0.5625	0.09375 0.625	0.125 0.5625	0.125 0.625	0.125 0.625	0.125 0.5625	0.15625 0.5625	0.125 0.625	0.15625 0.5625	0.15625 0.625	0.15625 0.625	0.15625 0.5625	0.1875 0.5625	0.15625 0.625	0.1875 0.5625	0.1875 0.625	0.1875 0.625	0.1875 0.5625	0.21875 0.5625	0.1875 0.625	0.21875 0.5625	0.21875 0.625	0.21875 0.625	0.21875 0.5625	0.25 0.5625	0.21875 0.625	0.25 0.5625	0.25 0.625	0.25 0.625	0.25 0.5625	0.28125 0.5625	0.25 0.625	0.28125 0.5625	0.28125 0.625	0.28125 0.625	0.28125 0.5625	0.3125 0.5625	0.28125 0.625	0.3125 0.5625	0.3125 0.625	0.3125 0.625	0.3125 0.5625	0.34375 0.5625	0.3125 0.625	0.34375 0.5625	0.34375 0.625	0.34375 0.625	0.34375 0.5625	0.375 0.5625	0.34375 0.625	0.375 0.5625	0.375 0.625	0.375 0.625	0.375 0.5625	0.40625 0.5625	0.375 0.625	0.40625 0.5625	0.40625 0.625	0.40625 0.625	0.40625 0.5625	0.4375 0.5625	0.40625 0.625	0.4375 0.5625	0.4375 0.625	0.4375 0.625	0.4375 0.5625	0.46875 0.5625	0.4375 0.625	0.46875 0.5625	0.46875 0.625	0.46875 0.625	0.46875 0.5625	0.5 0.5625	0.46875 0.625	0.5 0.5625	0.5 0.625	0.5 0.625	0.5 0.5625	0.53125 0.5625	0.5 0.625	0.53125 0.5625	0.53125 0.625	0.53125 0.625	0.53125 0.5625	0.5625 0.5625	0.53125 0.625	0.5625 0.5625	0.5625 0.625	0.5625 0.625	0.5625 0.5625	0.59375 0.5625	0.5625 0.625	0.59375 0.5625	0.59375 0.625	0.59375 0.625	0.59375 0.5625	0.625 0.5625	0.59375 0.625	0.625 0.5625	0.625 0.625	0.625 0.625	0.625 0.5625	0.65625 0.5625	0.625 0.625	0.65625 0.5625	0.65625 0.625	0.65625 0.625	0.65625 0.5625	0.6875 0.5625	0.65625 0.625	0.6875 0.5625	0.6875 0.625	0.6875 0.625	0.6875 0.5625	0.71875 0.5625	0.6875 0.625	0.71875 0.5625	0.71875 0.625	0.71875 0.625	0.71875 0.5625	0.75 0.5625	0.71875 0.625	0.75 0.5625	0.75 0.625	0.75 0.625	0.75 0.5625	0.78125 0.5625	0.75 0.625	0.78125 0.5625	0.78125 0.625	0.78125 0.625	0.78125 0.5625	0.8125 0.5625	0.78125 0.625	0.8125 0.5625	0.8125 0.625	0.8125 0.625	0.8125 0.5625	0.84375 0.5625	0.8125 0.625	0.84375 0.5625	0.84375 0.625	0.84375 0.625	0.84375 0.5625	0.875 0.5625	0.84375 0.625	0.875 0.5625	0.875 0.625	0.875 0.625	0.875 0.5625	0.90625 0.5625	0.875 0.625	0.90625 0.5625	0.90625 0.625	0.90625 0.625	0.90625 0.5625	0.9375 0.5625	0.90625 0.625	0.9375 0.5625	0.9375 0.625	0.9375 0.625	0.9375 0.5625	0.96875 0.5625	0.9375 0.625	0.96875 0.5625	0.96875 0.625	0.96875 0.625	0.96875 0.5625	1.0 0.5625	0.96875 0.625	1.0 0.5625	1.0 0.625	0.0 0.5625	0.0 0.5	0.03125 0.5	0.0 0.5625	0.03125 0.5	0.03125 0.5625	0.03125 0.5625	0.03125 0.5	0.0625 0.5	0.03125 0.5625	0.0625 0.5	0.0625 0.5625	0.0625 0.5625	0.0625 0.5	0.09375 0.5	0.0625 0.5625	0.09375 0.5	0.09375 0.5625	0.09375 0.5625	0.09375 0.5	0.125 0.5	0.09375 0.5625	0.125 0.5	0.125 0.5625	0.125 0.5625	0.125 0.5	0.15625 0.5	0.125 0.5625	0.15625 0.5	0.15625 0.5625	0.15625 0.5625	0.15625 0.5	0.1875 0.5	0.15625 0.5625	0.1875 0.5	0.1875 0.5625	0.1875 0.5625	0.1875 0.5	0.21875 0.5	0.1875 0.5625	0.21875 0.5	0.21875 0.5625	0.21875 0.5625	0.21875 0.5	0.25 0.5	0.21875 0.5625	0.25 0.5	0.25 0.5625	0.25 0.5625	0.25 0.5	0.28125 0.5	0.25 0.5625	0.28125 0.5	0.28125 0.5625	0.28125 0.5625	0.28125 0.5	0.3125 0.5	0.28125 0.5625	0.3125 0.5	0.3125 0.5625	0.3125 0.5625	0.3125 0.5	0.34375 0.5	0.3125 0.5625	0.34375 0.5	0.34375 0.5625	0.34375 0.5625	0.34375 0.5	0.375 0.5	0.34375 0.5625	0.375 0.5	0.375 0.5625	0.375 0.5625	0.375 0.5	0.40625 0.5	0.375 0.5625	0.40625 0.5	0.40625 0.5625	0.40625 0.5625	0.40625 0.5	0.4375 0.5	0.40625 0.5625	0.4375 0.5	0.4375 0.5625	0.4375 0.5625	0.4375 0.5	0.46875 0.5	0.4375 0.5625	0.46875 0.5	0.46875 0.5625	0.46875 0.5625	0.46875 0.5	0.5 0.5	0.46875 0.5625	0.5 0.5	0.5 0.5625	0.5 0.5625	0.5 0.5	0.53125 0.5	0.5 0.5625	0.53125 0.5	0.53125 0.5625	0.53125 0.5625	0.53125 0.5	0.5625 0.5	0.53125 0.5625	0.5625 0.5	0.5625 0.5625	0.5625 0.5625	0.5625 0.5	0.59375 0.5	0.5625 0.5625	0.59375 0.5	0.59375 0.5625	0.59375 0.5625	0.59375 0.5	0.625 0.5	0.59375 0.5625	0.625 0.5	0.625 0.5625	0.625 0.5625	0.625 0.5	0.65625 0.5	0.625 0.5625	0.65625 0.5	0.65625 0.5625	0.65625 0.5625	0.65625 0.5	0.6875 0.5	0.65625 0.5625	0.6875 0.5	0.6875 0.5625	0.6875 0.5625	0.6875 0.5	0.71875 0.5	0.6875 0.5625	0.71875 0.5	0.71875 0.5625	0.71875 0.5625	0.71875 0.5	0.75 0.5	0.71875 0.5625	0.75 0.5	0.75 0.5625	0.75 0.5625	0.75 0.5	0.78125 0.5	0.75 0.5625	0.78125 0.5	0.78125 0.5625	0.78125 0.5625	0.78125 0.5	0.8125 0.5	0.78125 0.5625	0.8125 0.5	0.8125 0.5625	0.8125 0.5625	0.8125 0.5	0.84375 0.5	0.8125 0.5625	0.84375 0.5	0.84375 0.5625	0.84375 0.5625	0.84375 0.5	0.875 0.5	0.84375 0.5625	0.875 0.5	0.875 0.5625	0.875 0.5625	0.875 0.5	0.90625 0.5	0.875 0.5625	0.90625 0.5	0.90625 0.5625	0.90625 0.5625	0.90625 0.5	0.9375 0.5	0.90625 0.5625	0.9375 0.5	0.9375 0.5625	0.9375 0.5625	0.9375 0.5	0.96875 0.5	0.9375 0.5625	0.96875 0.5	0.96875 0.5625	0.96875 0.5625	0.96875 0.5	1.0 0.5	0.96875 0.5625	1.0 0.5	1.0 0.5625	0.0 0.5	0.0 0.4375	0.03125 0.4375	0.0 0.5	0.03125 0.4375	0.03125 0.5	0.03125 0.5	0.03125 0.4375	0.0625 0.4375	0.03125 0.5	0.0625 0.4375	0.0625 0.5	0.0625 0.5	0.0625 0.4375	0.09375 0.4375	0.0625 0.5	0.09375 0.4375	0.09375 0.5	0.09375 0.5	0.09375 0.4375	0.125 0.4375	0.09375 0.5	0.125 0.4375	0.125 0.5	0.125 0.5	0.125 0.4375	0.15625 0.4375	0.125 0.5	0.15625 0.4375	0.15625 0.5	0.15625 0.5	0.15625 0.4375	0.1875 0.4375	0.15625 0.5	0.1875 0.4375	0.1875 0.5	0.1875 0.5	0.1875 0.4375	0.21875 0.4375	0.1875 0.5	0.21875 0.4375	0.21875 0.5	0.21875 0.5	0.21875 0.4375	0.25 0.4375	0.21875 0.5	0.25 0.4375	0.25 0.5	0.25 0.5	0.25 0.4375	0.28125 0.4375	0.25 0.5	0.28125 0.4375	0.28125 0.5	0.28125 0.5	0.28125 0.4375	0.3125 0.4375	0.28125 0.5	0.3125 0.4375	0.3125 0.5	0.3125 0.5	0.3125 0.4375	0.34375 0.4375	0.3125 0.5	0.34375 0.4375	0.34375 0.5	0.34375 0.5	0.34375 0.4375	0.375 0.4375	0.34375 0.5	0.375 0.4375	0.375 0.5	0.375 0.5	0.375 0.4375	0.40625 0.4375	0.375 0.5	0.40625 0.4375	0.40625 0.5	0.40625 0.5	0.40625 0.4375	0.4375 0.4375	0.40625 0.5	0.4375 0.4375	0.4375 0.5	0.4375 0.5	0.4375 0.4375	0.46875 0.4375	0.4375 0.5	0.46875 0.4375	0.46875 0.5	0.46875 0.5	0.46875 0.4375	0.5 0.4375	0.46875 0.5	0.5 0.4375	0.5 0.5	0.5 0.5	0.5 0.4375	0.53125 0.4375	0.5 0.5	0.53125 0.4375	0.53125 0.5	0.53125 0.5	0.53125 0.4375	0.5625 0.4375	0.53125 0.5	0.5625 0.4375	0.5625 0.5	0.5625 0.5	0.5625 0.4375	0.59375 0.4375	0.5625 0.5	0.59375 0.4375	0.59375 0.5	0.59375 0.5	0.59375 0.4375	0.625 0.4375	0.59375 0.5	0.625 0.4375	0.625 0.5	0.625 0.5	0.625 0.4375	0.65625 0.4375	0.625 0.5	0.65625 0.4375	0.65625 0.5	0.65625 0.5	0.65625 0.4375	0.6875 0.4375	0.65625 0.5	0.6875 0.4375	0.6875 0.5	0.6875 0.5	0.6875 0.4375	0.71875 0.4375	0.6875 0.5	0.71875 0.4375	0.71875 0.5	0.71875 0.5	0.71875 0.4375	0.75 0.4375	0.71875 0.5	0.75 0.4375	0.75 0.5	0.75 0.5	0.75 0.4375	0.78125 0.4375	0.75 0.5	0.78125 0.4375	0.78125 0.5	0.78125 0.5	0.78125 0.4375	0.8125 0.4375	0.78125 0.5	0.8125 0.4375	0.8125 0.5	0.8125 0.5	0.8125 0.4375	0.84375 0.4375	0.8125 0.5	0.84375 0.4375	0.84375 0.5	0.84375 0.5	0.84375 0.4375	0.875 0.4375	0.84375 0.5	0.875 0.4375	0.875 0.5	0.875 0.5	0.875 0.4375	0.90625 0.4375	0.875 0.5	0.90625 0.4375	0.90625 0.5	0.90625 0.5	0.90625 0.4375	0.9375 0.4375	0.90625 0.5	0.9375 0.4375	0.9375 0.5	0.9375 0.5	0.9375 0.4375	0.96875 0.4375	0.9375 0.5	0.96875 0.4375	0.96875 0.5	0.96875 0.5	0.96875 0.4375	1.0 0.4375	0.96875 0.5	1.0 0.4375	1.0 0.5	0.0 0.4375	0.0 0.375	0.03125 0.375	0.0 0.4375	0.03125 0.375	0.03125 0.4375	0.03125 0.4375	0.03125 0.375	0.0625 0.375	0.03125 0.4375	0.0625 0.375	0.0625 0.4375	0.0625 0.4375	0.0625 0.375	0.09375 0.375	0.0625 0.4375	0.09375 0.375	0.09375 0.4375	0.09375 0.4375	0.09375 0.375	0.125 0.375	0.09375 0.4375	0.125 0.375	0.125 0.4375	0.125 0.4375	0.125 0.375	0.15625 0.375	0.125 0.4375	0.15625 0.375	0.15625 0.4375	0.15625 0.4375	0.15625 0.375	0.1875 0.375	0.15625 0.4375	0.1875 0.375	0.1875 0.4375	0.1875 0.4375	0.1875 0.375	0.21875 0.375	0.1875 0.4375	0.21875 0.375	0.21875 0.4375	0.21875 0.4375	0.21875 0.375	0.25 0.375	0.21875 0.4375	0.25 0.375	0.25 0.4375	0.25 0.4375	0.25 0.375	0.28125 0.375	0.25 0.4375	0.28125 0.375	0.28125 0.4375	0.28125 0.4375	0.28125 0.375	0.3125 0.375	0.28125 0.4375	0.3125 0.375	0.3125 0.4375	0.3125 0.4375	0.3125 0.375	0.34375 0.375	0.3125 0.4375	0.34375 0.375	0.34375 0.4375	0.34375 0.4375	0.34375 0.375	0.375 0.375	0.34375 0.4375	0.375 0.375	0.375 0.4375	0.375 0.4375	0.375 0.375	0.40625 0.375	0.375 0.4375	0.40625 0.375	0.40625 0.4375	0.40625 0.4375	0.40625 0.375	0.4375 0.375	0.40625 0.4375	0.4375 0.375	0.4375 0.4375	0.4375 0.4375	0.4375 0.375	0.46875 0.375	0.4375 0.4375	0.46875 0.375	0.46875 0.4375	0.46875 0.4375	0.46875 0.375	0.5 0.375	0.46875 0.4375	0.5 0.375	0.5 0.4375	0.5 0.4375	0.5 0.375	0.53125 0.375	0.5 0.4375	0.53125 0.375	0.53125 0.4375	0.53125 0.4375	0.53125 0.375	0.5625 0.375	0.53125 0.4375	0.5625 0.375	0.5625 0.4375	0.5625 0.4375	0.5625 0.375	0.59375 0.375	0.5625 0.4375	0.59375 0.375	0.59375 0.4375	0.59375 0.4375	0.59375 0.375	0.625 0.375	0.59375 0.4375	0.625 0.375	0.625 0.4375	0.625 0.4375	0.625 0.375	0.65625 0.375	0.625 0.4375	0.65625 0.375	0.65625 0.4375	0.65625 0.4375	0.65625 0.375	0.6875 0.375	0.65625 0.4375	0.6875 0.375	0.6875 0.4375	0.6875 0.4375	0.6875 0.375	0.71875 0.375	0.6875 0.4375	0.71875 0.375	0.71875 0.4375	0.71875 0.4375	0.71875 0.375	0.75 0.375	0.71875 0.4375	0.75 0.375	0.75 0.4375	0.75 0.4375	0.75 0.375	0.78125 0.375	0.75 0.4375	0.78125 0.375	0.78125 0.4375	0.78125 0.4375	0.78125 0.375	0.8125 0.375	0.78125 0.4375	0.8125 0.375	0.8125 0.4375	0.8125 0.4375	0.8125 0.375	0.84375 0.375	0.8125 0.4375	0.84375 0.375	0.84375 0.4375	0.84375 0.4375	0.84375 0.375	0.875 0.375	0.84375 0.4375	0.875 0.375	0.875 0.4375	0.875 0.4375	0.875 0.375	0.90625 0.375	0.875 0.4375	0.90625 0.375	0.90625 0.4375	0.90625 0.4375	0.90625 0.375	0.9375 0.375	0.90625 0.4375	0.9375 0.375	0.9375 0.4375	0.9375 0.4375	0.9375 0.375	0.96875 0.375	0.9375 0.4375	0.96875 0.375	0.96875 0.4375	0.96875 0.4375	0.96875 0.375	1.0 0.375	0.96875 0.4375	1.0 0.375	1.0 0.4375	0.0 0.375	0.0 0.3125	0.03125 0.3125	0.0 0.375	0.03125 0.3125	0.03125 0.375	0.03125 0.375	0.03125 0.3125	0.0625 0.3125	0.03125 0.375	0.0625 0.3125	0.0625 0.375	0.0625 0.375	0.0625 0.3125	0.09375 0.3125	0.0625 0.375	0.09375 0.3125	0.09375 0.375	0.09375 0.375	0.09375 0.3125	0.125 0.3125	0.09375 0.375	0.125 0.3125	0.125 0.375	0.125 0.375	0.125 0.3125	0.15625 0.3125	0.125 0.375	0.15625 0.3125	0.15625 0.375	0.15625 0.375	0.15625 0.3125	0.1875 0.3125	0.15625 0.375	0.1875 0.3125	0.1875 0.375	0.1875 0.375	0.1875 0.3125	0.21875 0.3125	0.1875 0.375	0.21875 0.3125	0.21875 0.375	0.21875 0.375	0.21875 0.3125	0.25 0.3125	0.21875 0.375	0.25 0.3125	0.25 0.375	0.25 0.375	0.25 0.3125	0.28125 0.3125	0.25 0.375	0.28125 0.3125	0.28125 0.375	0.28125 0.375	0.28125 0.3125	0.3125 0.3125	0.28125 0.375	0.3125 0.3125	0.3125 0.375	0.3125 0.375	0.3125 0.3125	0.34375 0.3125	0.3125 0.375	0.34375 0.3125	0.34375 0.375	0.34375 0.375	0.34375 0.3125	0.375 0.3125	0.34375 0.375	0.375 0.3125	0.375 0.375	0.375 0.375	0.375 0.3125	0.40625 0.3125	0.375 0.375	0.40625 0.3125	0.40625 0.375	0.40625 0.375	0.40625 0.3125	0.4375 0.3125	0.40625 0.375	0.4375 0.3125	0.4375 0.375	0.4375 0.375	0.4375 0.3125	0.46875 0.3125	0.4375 0.375	0.46875 0.3125	0.46875 0.375	0.46875 0.375	0.46875 0.3125	0.5 0.3125	0.46875 0.375	0.5 0.3125	0.5 0.375	0.5 0.375	0.5 0.3125	0.53125 0.3125	0.5 0.375	0.53125 0.3125	0.53125 0.375	0.53125 0.375	0.53125 0.3125	0.5625 0.3125	0.53125 0.375	0.5625 0.3125	0.5625 0.375	0.5625 0.375	0.5625 0.3125	0.59375 0.3125	0.5625 0.375	0.59375 0.3125	0.59375 0.375	0.59375 0.375	0.59375 0.3125	0.625 0.3125	0.59375 0.375	0.625 0.3125	0.625 0.375	0.625 0.375	0.625 0.3125	0.65625 0.3125	0.625 0.375	0.65625 0.3125	0.65625 0.375	0.65625 0.375	0.65625 0.3125	0.6875 0.3125	0.65625 0.375	0.6875 0.3125	0.6875 0.375	0.6875 0.375	0.6875 0.3125	0.71875 0.3125	0.6875 0.375	0.71875 0.3125	0.71875 0.375	0.71875 0.375	0.71875 0.3125	0.75 0.3125	0.71875 0.375	0.75 0.3125	0.75 0.375	0.75 0.375	0.75 0.3125	0.78125 0.3125	0.75 0.375	0.78125 0.3125	0.78125 0.375	0.78125 0.375	0.78125 0.3125	0.8125 0.3125	0.78125 0.375	0.8125 0.3125	0.8125 0.375	0.8125 0.375	0.8125 0.3125	0.84375 0.3125	0.8125 0.375	0.84375 0.3125	0.84375 0.375	0.84375 0.375	0.84375 0.3125	0.875 0.3125	0.84375 0.375	0.875 0.3125	0.875 0.375	0.875 0.375	0.875 0.3125	0.90625 0.3125	0.875 0.375	0.90625 0.3125	0.90625 0.375	0.90625 0.375	0.90625 0.3125	0.9375 0.3125	0.90625 0.375	0.9375 0.3125	0.9375 0.375	0.9375 0.375	0.9375 0.3125	0.96875 0.3125	0.9375 0.375	0.96875 0.3125	0.96875 0.375	0.96875 0.375	0.96875 0.3125	1.0 0.3125	0.96875 0.375	1.0 0.3125	1.0 0.375	0.0 0.3125	0.0 0.25	0.03125 0.25	0.0 0.3125	0.03125 0.25	0.03125 0.3125	0.03125 0.3125	0.03125 0.25	0.0625 0.25	0.03125 0.3125	0.0625 0.25	0.0625 0.3125	0.0625 0.3125	0.0625 0.25	0.09375 0.25	0.0625 0.3125	0.09375 0.25	0.09375 0.3125	0.09375 0.3125	0.09375 0.25	0.125 0.25	0.09375 0.3125	0.125 0.25	0.125 0.3125	0.125 0.3125	0.125 0.25	0.15625 0.25	0.125 0.3125	0.15625 0.25	0.15625 0.3125	0.15625 0.3125	0.15625 0.25	0.1875 0.25	0.15625 0.3125	0.1875 0.25	0.1875 0.3125	0.1875 0.3125	0.1875 0.25	0.21875 0.25	0.1875 0.3125	0.21875 0.25	0.21875 0.3125	0.21875 0.3125	0.21875 0.25	0.25 0.25	0.21875 0.3125	0.25 0.25	0.25 0.3125	0.25 0.3125	0.25 0.25	0.28125 0.25	0.25 0.3125	0.28125 0.25	0.28125 0.3125	0.28125 0.3125	0.28125 0.25	0.3125 0.25	0.28125 0.3125	0.3125 0.25	0.3125 0.3125	0.3125 0.3125	0.3125 0.25	0.34375 0.25	0.3125 0.3125	0.34375 0.25	0.34375 0.3125	0.34375 0.3125	0.34375 0.25	0.375 0.25	0.34375 0.3125	0.375 0.25	0.375 0.3125	0.375 0.3125	0.375 0.25	0.40625 0.25	0.375 0.3125	0.40625 0.25	0.40625 0.3125	0.40625 0.3125	0.40625 0.25	0.4375 0.25	0.40625 0.3125	0.4375 0.25	0.4375 0.3125	0.4375 0.3125	0.4375 0.25	0.46875 0.25	0.4375 0.3125	0.46875 0.25	0.46875 0.3125	0.46875 0.3125	0.46875 0.25	0.5 0.25	0.46875 0.3125	0.5 0.25	0.5 0.3125	0.5 0.3125	0.5 0.25	0.53125 0.25	0.5 0.3125	0.53125 0.25	0.53125 0.3125	0.53125 0.3125	0.53125 0.25	0.5625 0.25	0.53125 0.3125	0.5625 0.25	0.5625 0.3125	0.5625 0.3125	0.5625 0.25	0.59375 0.25	0.5625 0.3125	0.59375 0.25	0.59375 0.3125	0.59375 0.3125	0.59375 0.25	0.625 0.25	0.59375 0.3125	0.625 0.25	0.625 0.3125	0.625 0.3125	0.625 0.25	0.65625 0.25	0.625 0.3125	0.65625 0.25	0.65625 0.3125	0.65625 0.3125	0.65625 0.25	0.6875 0.25	0.65625 0.3125	0.6875 0.25	0.6875 0.3125	0.6875 0.3125	0.6875 0.25	0.71875 0.25	0.6875 0.3125	0.71875 0.25	0.71875 0.3125	0.71875 0.3125	0.71875 0.25	0.75 0.25	0.71875 0.3125	0.75 0.25	0.75 0.3125	0.75 0.3125	0.75 0.25	0.78125 0.25	0.75 0.3125	0.78125 0.25	0.78125 0.3125	0.78125 0.3125	0.78125 0.25	0.8125 0.25	0.78125 0.3125	0.8125 0.25	0.8125 0.3125	0.8125 0.3125	0.8125 0.25	0.84375 0.25	0.8125 0.3125	0.84375 0.25	0.84375 0.3125	0.84375 0.3125	0.84375 0.25	0.875 0.25	0.84375 0.3125	0.875 0.25	0.875 0.3125	0.875 0.3125	0.875 0.25	0.90625 0.25	0.875 0.3125	0.90625 0.25	0.90625 0.3125	0.90625 0.3125	0.90625 0.25	0.9375 0.25	0.90625 0.3125	0.9375 0.25	0.9375 0.3125	0.9375 0.3125	0.9375 0.25	0.96875 0.25	0.9375 0.3125	0.96875 0.25	0.96875 0.3125	0.96875 0.3125	0.96875 0.25	1.0 0.25	0.96875 0.3125	1.0 0.25	1.0 0.3125	0.0 0.25	0.0 0.1875	0.03125 0.1875	0.0 0.25	0.03125 0.1875	0.03125 0.25	0.03125 0.25	0.03125 0.1875	0.0625 0.1875	0.03125 0.25	0.0625 0.1875	0.0625 0.25	0.0625 0.25	0.0625 0.1875	0.09375 0.1875	0.0625 0.25	0.09375 0.1875	0.09375 0.25	0.09375 0.25	0.09375 0.1875	0.125 0.1875	0.09375 0.25	0.125 0.1875	0.125 0.25	0.125 0.25	0.125 0.1875	0.15625 0.1875	0.125 0.25	0.15625 0.1875	0.15625 0.25	0.15625 0.25	0.15625 0.1875	0.1875 0.1875	0.15625 0.25	0.1875 0.1875	0.1875 0.25	0.1875 0.25	0.1875 0.1875	0.21875 0.1875	0.1875 0.25	0.21875 0.1875	0.21875 0.25	0.21875 0.25	0.21875 0.1875	0.25 0.1875	0.21875 0.25	0.25 0.1875	0.25 0.25	0.25 0.25	0.25 0.1875	0.28125 0.1875	0.25 0.25	0.28125 0.1875	0.28125 0.25	0.28125 0.25	0.28125 0.1875	0.3125 0.1875	0.28125 0.25	0.3125 0.1875	0.3125 0.25	0.3125 0.25	0.3125 0.1875	0.34375 0.1875	0.3125 0.25	0.34375 0.1875	0.34375 0.25	0.34375 0.25	0.34375 0.1875	0.375 0.1875	0.34375 0.25	0.375 0.1875	0.375 0.25	0.375 0.25	0.375 0.1875	0.40625 0.1875	0.375 0.25	0.40625 0.1875	0.40625 0.25	0.40625 0.25	0.40625 0.1875	0.4375 0.1875	0.40625 0.25	0.4375 0.1875	0.4375 0.25	0.4375 0.25	0.4375 0.1875	0.46875 0.1875	0.4375 0.25	0.46875 0.1875	0.46875 0.25	0.46875 0.25	0.46875 0.1875	0.5 0.1875	0.46875 0.25	0.5 0.1875	0.5 0.25	0.5 0.25	0.5 0.1875	0.53125 0.1875	0.5 0.25	0.53125 0.1875	0.53125 0.25	0.53125 0.25	0.53125 0.1875	0.5625 0.1875	0.53125 0.25	0.5625 0.1875	0.5625 0.25	0.5625 0.25	0.5625 0.1875	0.59375 0.1875	0.5625 0.25	0.59375 0.1875	0.59375 0.25	0.59375 0.25	0.59375 0.1875	0.625 0.1875	0.59375 0.25	0.625 0.1875	0.625 0.25	0.625 0.25	0.625 0.1875	0.65625 0.1875	0.625 0.25	0.65625 0.1875	0.65625 0.25	0.65625 0.25	0.65625 0.1875	0.6875 0.1875	0.65625 0.25	0.6875 0.1875	0.6875 0.25	0.6875 0.25	0.6875 0.1875	0.71875 0.1875	0.6875 0.25	0.71875 0.1875	0.71875 0.25	0.71875 0.25	0.71875 0.1875	0.75 0.1875	0.71875 0.25	0.75 0.1875	0.75 0.25	0.75 0.25	0.75 0.1875	0.78125 0.1875	0.75 0.25	0.78125 0.1875	0.78125 0.25	0.78125 0.25	0.78125 0.1875	0.8125 0.1875	0.78125 0.25	0.8125 0.1875	0.8125 0.25	0.8125 0.25	0.8125 0.1875	0.84375 0.1875	0.8125 0.25	0.84375 0.1875	0.84375 0.25	0.84375 0.25	0.84375 0.1875	0.875 0.1875	0.84375 0.25	0.875 0.1875	0.875 0.25	0.875 0.25	0.875 0.1875	0.90625 0.1875	0.875 0.25	0.90625 0.1875	0.90625 0.25	0.90625 0.25	0.90625 0.1875	0.9375 0.1875	0.90625 0.25	0.9375 0.1875	0.9375 0.25	0.9375 0.25	0.9375 0.1875	0.96875 0.1875	0.9375 0.25	0.96875 0.1875	0.96875 0.25	0.96875 0.25	0.96875 0.1875	1.0 0.1875	0.96875 0.25	1.0 0.1875	1.0 0.25	0.0 0.1875	0.0 0.125	0.03125 0.125	0.0 0.1875	0.03125 0.125	0.03125 0.1875	0.03125 0.1875	0.03125 0.125	0.0625 0.125	0.03125 0.1875	0.0625 0.125	0.0625 0.1875	0.0625 0.1875	0.0625 0.125	0.09375 0.125	0.0625 0.1875	0.09375 0.125	0.09375 0.1875	0.09375 0.1875	0.09375 0.125	0.125 0.125	0.09375 0.1875	0.125 0.125	0.125 0.1875	0.125 0.1875	0.125 0.125	0.15625 0.125	0.125 0.1875	0.15625 0.125	0.15625 0.1875	0.15625 0.1875	0.15625 0.125	0.1875 0.125	0.15625 0.1875	0.1875 0.125	0.1875 0.1875	0.1875 0.1875	0.1875 0.125	0.21875 0.125	0.1875 0.1875	0.21875 0.125	0.21875 0.1875	0.21875 0.1875	0.21875 0.125	0.25 0.125	0.21875 0.1875	0.25 0.125	0.25 0.1875	0.25 0.1875	0.25 0.125	0.28125 0.125	0.25 0.1875	0.28125 0.125	0.28125 0.1875	0.28125 0.1875	0.28125 0.125	0.3125 0.125	0.28125 0.1875	0.3125 0.125	0.3125 0.1875	0.3125 0.1875	0.3125 0.125	0.34375 0.125	0.3125 0.1875	0.34375 0.125	0.34375 0.1875	0.34375 0.1875	0.34375 0.125	0.375 0.125	0.34375 0.1875	0.375 0.125	0.375 0.1875	0.375 0.1875	0.375 0.125	0.40625 0.125	0.375 0.1875	0.40625 0.125	0.40625 0.1875	0.40625 0.1875	0.40625 0.125	0.4375 0.125	0.40625 0.1875	0.4375 0.125	0.4375 0.1875	0.4375 0.1875	0.4375 0.125	0.46875 0.125	0.4375 0.1875	0.46875 0.125	0.46875 0.1875	0.46875 0.1875	0.46875 0.125	0.5 0.125	0.46875 0.1875	0.5 0.125	0.5 0.1875	0.5 0.1875	0.5 0.125	0.53125 0.125	0.5 0.1875	0.53125 0.125	0.53125 0.1875	0.53125 0.1875	0.53125 0.125	0.5625 0.125	0.53125 0.1875	0.5625 0.125	0.5625 0.1875	0.5625 0.1875	0.5625 0.125	0.59375 0.125	0.5625 0.1875	0.59375 0.125	0.59375 0.1875	0.59375 0.1875	0.59375 0.125	0.625 0.125	0.59375 0.1875	0.625 0.125	0.625 0.1875	0.625 0.1875	0.625 0.125	0.65625 0.125	0.625 0.1875	0.65625 0.125	0.65625 0.1875	0.65625 0.1875	0.65625 0.125	0.6875 0.125	0.65625 0.1875	0.6875 0.125	0.6875 0.1875	0.6875 0.1875	0.6875 0.125	0.71875 0.125	0.6875 0.1875	0.71875 0.125	0.71875 0.1875	0.71875 0.1875	0.71875 0.125	0.75 0.125	0.71875 0.1875	0.75 0.125	0.75 0.1875	0.75 0.1875	0.75 0.125	0.78125 0.125	0.75 0.1875	0.78125 0.125	0.78125 0.1875	0.78125 0.1875	0.78125 0.125	0.8125 0.125	0.78125 0.1875	0.8125 0.125	0.8125 0.1875	0.8125 0.1875	0.8125 0.125	0.84375 0.125	0.8125 0.1875	0.84375 0.125	0.84375 0.1875	0.84375 0.1875	0.84375 0.125	0.875 0.125	0.84375 0.1875	0.875 0.125	0.875 0.1875	0.875 0.1875	0.875 0.125	0.90625 0.125	0.875 0.1875	0.90625 0.125	0.90625 0.1875	0.90625 0.1875	0.90625 0.125	0.9375 0.125	0.90625 0.1875	0.9375 0.125	0.9375 0.1875	0.9375 0.1875	0.9375 0.125	0.96875 0.125	0.9375 0.1875	0.96875 0.125	0.96875 0.1875	0.96875 0.1875	0.96875 0.125	1.0 0.125	0.96875 0.1875	1.0 0.125	1.0 0.1875	0.0 0.125	0.0 0.0624999	0.03125 0.0624999	0.0 0.125	0.03125 0.0624999	0.03125 0.125	0.03125 0.125	0.03125 0.0624999	0.0625 0.0624999	0.03125 0.125	0.0625 0.0624999	0.0625 0.125	0.0625 0.125	0.0625 0.0624999	0.09375 0.0624999	0.0625 0.125	0.09375 0.0624999	0.09375 0.125	0.09375 0.125	0.09375 0.0624999	0.125 0.0624999	0.09375 0.125	0.125 0.0624999	0.125 0.125	0.125 0.125	0.125 0.0624999	0.15625 0.0624999	0.125 0.125	0.15625 0.0624999	0.15625 0.125	0.15625 0.125	0.15625 0.0624999	0.1875 0.0624999	0.15625 0.125	0.1875 0.0624999	0.1875 0.125	0.1875 0.125	0.1875 0.0624999	0.21875 0.0624999	0.1875 0.125	0.21875 0.0624999	0.21875 0.125	0.21875 0.125	0.21875 0.0624999	0.25 0.0624999	0.21875 0.125	0.25 0.0624999	0.25 0.125	0.25 0.125	0.25 0.0624999	0.28125 0.0624999	0.25 0.125	0.28125 0.0624999	0.28125 0.125	0.28125 0.125	0.28125 0.0624999	0.3125 0.0624999	0.28125 0.125	0.3125 0.0624999	0.3125 0.125	0.3125 0.125	0.3125 0.0624999	0.34375 0.0624999	0.3125 0.125	0.34375 0.0624999	0.34375 0.125	0.34375 0.125	0.34375 0.0624999	0.375 0.0624999	0.34375 0.125	0.375 0.0624999	0.375 0.125	0.375 0.125	0.375 0.0624999	0.40625 0.0624999	0.375 0.125	0.40625 0.0624999	0.40625 0.125	0.40625 0.125	0.40625 0.0624999	0.4375 0.0624999	0.40625 0.125	0.4375 0.0624999	0.4375 0.125	0.4375 0.125	0.4375 0.0624999	0.46875 0.0624999	0.4375 0.125	0.46875 0.0624999	0.46875 0.125	0.46875 0.125	0.46875 0.0624999	0.5 0.0624999	0.46875 0.125	0.5 0.0624999	0.5 0.125	0.5 0.125	0.5 0.0624999	0.53125 0.0624999	0.5 0.125	0.53125 0.0624999	0.53125 0.125	0.53125 0.125	0.53125 0.0624999	0.5625 0.0624999	0.53125 0.125	0.5625 0.0624999	0.5625 0.125	0.5625 0.125	0.5625 0.0624999	0.59375 0.0624999	0.5625 0.125	0.59375 0.0624999	0.59375 0.125	0.59375 0.125	0.59375 0.0624999	0.625 0.0624999	0.59375 0.125	0.625 0.0624999	0.625 0.125	0.625 0.125	0.625 0.0624999	0.65625 0.0624999	0.625 0.125	0.65625 0.0624999	0.65625 0.125	0.65625 0.125	0.65625 0.0624999	0.6875 0.0624999	0.65625 0.125	0.6875 0.0624999	0.6875 0.125	0.6875 0.125	0.6875 0.0624999	0.71875 0.0624999	0.6875 0.125	0.71875 0.0624999	0.71875 0.125	0.71875 0.125	0.71875 0.0624999	0.75 0.0624999	0.71875 0.125	0.75 0.0624999	0.75 0.125	0.75 0.125	0.75 0.0624999	0.78125 0.0624999	0.75 0.125	0.78125 0.0624999	0.78125 0.125	0.78125 0.125	0.78125 0.0624999	0.8125 0.0624999	0.78125 0.125	0.8125 0.0624999	0.8125 0.125	0.8125 0.125	0.8125 0.0624999	0.84375 0.0624999	0.8125 0.125	0.84375 0.0624999	0.84375 0.125	0.84375 0.125	0.84375 0.0624999	0.875 0.0624999	0.84375 0.125	0.875 0.0624999	0.875 0.125	0.875 0.125	0.875 0.0624999	0.90625 0.0624999	0.875 0.125	0.90625 0.0624999	0.90625 0.125	0.90625 0.125	0.90625 0.0624999	0.9375 0.0624999	0.90625 0.125	0.9375 0.0624999	0.9375 0.125	0.9375 0.125	0.9375 0.0624999	0.96875 0.0624999	0.9375 0.125	0.96875 0.0624999	0.96875 0.125	0.96875 0.125	0.96875 0.0624999	1.0 0.0624999	0.96875 0.125	1.0 0.0624999	1.0 0.125	0.0 -1.19209e-007	0.03125 0.0624999	0.0 0.0624999	0.03125 -1.19209e-007	0.0625 0.0624999	0.03125 0.0624999	0.0625 -1.19209e-007	0.09375 0.0624999	0.0625 0.0624999	0.09375 -1.19209e-007	0.125 0.0624999	0.09375 0.0624999	0.125 -1.19209e-007	0.15625 0.0624999	0.125 0.0624999	0.15625 -1.19209e-007	0.1875 0.0624999	0.15625 0.0624999	0.1875 -1.19209e-007	0.21875 0.0624999	0.1875 0.0624999	0.21875 -1.19209e-007	0.25 0.0624999	0.21875 0.0624999	0.25 -1.19209e-007	0.28125 0.0624999	0.25 0.0624999	0.28125 -1.19209e-007	0.3125 0.0624999	0.28125 0.0624999	0.3125 -1.19209e-007	0.34375 0.0624999	0.3125 0.0624999	0.34375 -1.19209e-007	0.375 0.0624999	0.34375 0.0624999	0.375 -1.19209e-007	0.40625 0.0624999	0.375 0.0624999	0.40625 -1.19209e-007	0.4375 0.0624999	0.40625 0.0624999	0.4375 -1.19209e-007	0.46875 0.0624999	0.4375 0.0624999	0.46875 -1.19209e-007	0.5 0.0624999	0.46875 0.0624999	0.5 -1.19209e-007	0.53125 0.0624999	0.5 0.0624999	0.53125 -1.19209e-007	0.5625 0.0624999	0.53125 0.0624999	0.5625 -1.19209e-007	0.59375 0.0624999	0.5625 0.0624999	0.59375 -1.19209e-007	0.625 0.0624999	0.59375 0.0624999	0.625 -1.19209e-007	0.65625 0.0624999	0.625 0.0624999	0.65625 -1.19209e-007	0.6875 0.0624999	0.65625 0.0624999	0.6875 -1.19209e-007	0.71875 0.0624999	0.6875 0.0624999	0.71875 -1.19209e-007	0.75 0.0624999	0.71875 0.0624999	0.75 -1.19209e-007	0.78125 0.0624999	0.75 0.0624999	0.78125 -1.19209e-007	0.8125 0.0624999	0.78125 0.0624999	0.8125 -1.19209e-007	0.84375 0.0624999	0.8125 0.0624999	0.84375 -1.19209e-007	0.875 0.0624999	0.84375 0.0624999	0.875 -1.19209e-007	0.90625 0.0624999	0.875 0.0624999	0.90625 -1.19209e-007	0.9375 0.0624999	0.90625 0.0624999	0.9375 -1.19209e-007	0.96875 0.0624999	0.9375 0.0624999	0.96875 -1.19209e-007	1.0 0.0624999	0.96875 0.0624999]\n" to:out_file
		format "\"normal N\" [0.0 0.0 1.0	0.0 0.200969 0.979598	-0.039207 0.197107 0.979598	0.0 0.0 1.0	-0.039207 0.197107 0.979598	-0.0769074 0.185671 0.979598	0.0 0.0 1.0	-0.0769074 0.185671 0.979598	-0.111652 0.1671 0.979598	0.0 0.0 1.0	-0.111652 0.1671 0.979598	-0.142106 0.142107 0.979598	0.0 0.0 1.0	-0.142106 0.142107 0.979598	-0.1671 0.111652 0.979598	0.0 0.0 1.0	-0.1671 0.111652 0.979598	-0.185671 0.0769074 0.979598	0.0 0.0 1.0	-0.185671 0.0769074 0.979598	-0.197107 0.039207 0.979598	0.0 0.0 1.0	-0.197107 0.039207 0.979598	-0.200969 0.0 0.979598	0.0 0.0 1.0	-0.200969 0.0 0.979598	-0.197107 -0.0392072 0.979598	0.0 0.0 1.0	-0.197107 -0.0392072 0.979598	-0.185671 -0.0769075 0.979598	0.0 0.0 1.0	-0.185671 -0.0769075 0.979598	-0.167099 -0.111652 0.979598	0.0 0.0 1.0	-0.167099 -0.111652 0.979598	-0.142106 -0.142107 0.979598	0.0 0.0 1.0	-0.142106 -0.142107 0.979598	-0.111652 -0.1671 0.979598	0.0 0.0 1.0	-0.111652 -0.1671 0.979598	-0.0769073 -0.185671 0.979598	0.0 0.0 1.0	-0.0769073 -0.185671 0.979598	-0.0392069 -0.197107 0.979598	0.0 0.0 1.0	-0.0392069 -0.197107 0.979598	1.81967e-007 -0.200969 0.979598	0.0 0.0 1.0	1.81967e-007 -0.200969 0.979598	0.0392073 -0.197107 0.979598	0.0 0.0 1.0	0.0392073 -0.197107 0.979598	0.0769077 -0.185671 0.979598	0.0 0.0 1.0	0.0769077 -0.185671 0.979598	0.111653 -0.167099 0.979598	0.0 0.0 1.0	0.111653 -0.167099 0.979598	0.142107 -0.142106 0.979598	0.0 0.0 1.0	0.142107 -0.142106 0.979598	0.1671 -0.111652 0.979598	0.0 0.0 1.0	0.1671 -0.111652 0.979598	0.185671 -0.0769072 0.979598	0.0 0.0 1.0	0.185671 -0.0769072 0.979598	0.197107 -0.0392068 0.979598	0.0 0.0 1.0	0.197107 -0.0392068 0.979598	0.200969 3.22768e-007 0.979598	0.0 0.0 1.0	0.200969 3.22768e-007 0.979598	0.197107 0.0392074 0.979598	0.0 0.0 1.0	0.197107 0.0392074 0.979598	0.185671 0.0769077 0.979598	0.0 0.0 1.0	0.185671 0.0769077 0.979598	0.167099 0.111653 0.979598	0.0 0.0 1.0	0.167099 0.111653 0.979598	0.142106 0.142107 0.979598	0.0 0.0 1.0	0.142106 0.142107 0.979598	0.111652 0.1671 0.979598	0.0 0.0 1.0	0.111652 0.1671 0.979598	0.0769071 0.185671 0.979598	0.0 0.0 1.0	0.0769071 0.185671 0.979598	0.0392067 0.197107 0.979598	0.0 0.0 1.0	0.0392067 0.197107 0.979598	0.0 0.200969 0.979598	0.0 0.200969 0.979598	-2.84465e-007 0.387887 0.921707	-0.075673 0.380434 0.921707	0.0 0.200969 0.979598	-0.075673 0.380434 0.921707	-0.039207 0.197107 0.979598	-0.039207 0.197107 0.979598	-0.075673 0.380434 0.921707	-0.148438 0.358361 0.921707	-0.039207 0.197107 0.979598	-0.148438 0.358361 0.921707	-0.0769074 0.185671 0.979598	-0.0769074 0.185671 0.979598	-0.148438 0.358361 0.921707	-0.215498 0.322516 0.921707	-0.0769074 0.185671 0.979598	-0.215498 0.322516 0.921707	-0.111652 0.1671 0.979598	-0.111652 0.1671 0.979598	-0.215498 0.322516 0.921707	-0.274277 0.274277 0.921707	-0.111652 0.1671 0.979598	-0.274277 0.274277 0.921707	-0.142106 0.142107 0.979598	-0.142106 0.142107 0.979598	-0.274277 0.274277 0.921707	-0.322516 0.215498 0.921707	-0.142106 0.142107 0.979598	-0.322516 0.215498 0.921707	-0.1671 0.111652 0.979598	-0.1671 0.111652 0.979598	-0.322516 0.215498 0.921707	-0.358361 0.148438 0.921707	-0.1671 0.111652 0.979598	-0.358361 0.148438 0.921707	-0.185671 0.0769074 0.979598	-0.185671 0.0769074 0.979598	-0.358361 0.148438 0.921707	-0.380434 0.0756729 0.921707	-0.185671 0.0769074 0.979598	-0.380434 0.0756729 0.921707	-0.197107 0.039207 0.979598	-0.197107 0.039207 0.979598	-0.380434 0.0756729 0.921707	-0.387887 0.0 0.921707	-0.197107 0.039207 0.979598	-0.387887 0.0 0.921707	-0.200969 0.0 0.979598	-0.200969 0.0 0.979598	-0.387887 0.0 0.921707	-0.380434 -0.0756732 0.921707	-0.200969 0.0 0.979598	-0.380434 -0.0756732 0.921707	-0.197107 -0.0392072 0.979598	-0.197107 -0.0392072 0.979598	-0.380434 -0.0756732 0.921707	-0.358361 -0.148438 0.921707	-0.197107 -0.0392072 0.979598	-0.358361 -0.148438 0.921707	-0.185671 -0.0769075 0.979598	-0.185671 -0.0769075 0.979598	-0.358361 -0.148438 0.921707	-0.322516 -0.215499 0.921707	-0.185671 -0.0769075 0.979598	-0.322516 -0.215499 0.921707	-0.167099 -0.111652 0.979598	-0.167099 -0.111652 0.979598	-0.322516 -0.215499 0.921707	-0.274277 -0.274278 0.921707	-0.167099 -0.111652 0.979598	-0.274277 -0.274278 0.921707	-0.142106 -0.142107 0.979598	-0.142106 -0.142107 0.979598	-0.274277 -0.274278 0.921707	-0.215498 -0.322516 0.921707	-0.142106 -0.142107 0.979598	-0.215498 -0.322516 0.921707	-0.111652 -0.1671 0.979598	-0.111652 -0.1671 0.979598	-0.215498 -0.322516 0.921707	-0.148438 -0.358361 0.921707	-0.111652 -0.1671 0.979598	-0.148438 -0.358361 0.921707	-0.0769073 -0.185671 0.979598	-0.0769073 -0.185671 0.979598	-0.148438 -0.358361 0.921707	-0.0756727 -0.380434 0.921707	-0.0769073 -0.185671 0.979598	-0.0756727 -0.380434 0.921707	-0.0392069 -0.197107 0.979598	-0.0392069 -0.197107 0.979598	-0.0756727 -0.380434 0.921707	3.96817e-007 -0.387887 0.921707	-0.0392069 -0.197107 0.979598	3.96817e-007 -0.387887 0.921707	1.81967e-007 -0.200969 0.979598	1.81967e-007 -0.200969 0.979598	3.96817e-007 -0.387887 0.921707	0.0756734 -0.380434 0.921707	1.81967e-007 -0.200969 0.979598	0.0756734 -0.380434 0.921707	0.0392073 -0.197107 0.979598	0.0392073 -0.197107 0.979598	0.0756734 -0.380434 0.921707	0.148438 -0.358361 0.921707	0.0392073 -0.197107 0.979598	0.148438 -0.358361 0.921707	0.0769077 -0.185671 0.979598	0.0769077 -0.185671 0.979598	0.148438 -0.358361 0.921707	0.215499 -0.322516 0.921707	0.0769077 -0.185671 0.979598	0.215499 -0.322516 0.921707	0.111653 -0.167099 0.979598	0.111653 -0.167099 0.979598	0.215499 -0.322516 0.921707	0.274278 -0.274277 0.921707	0.111653 -0.167099 0.979598	0.274278 -0.274277 0.921707	0.142107 -0.142106 0.979598	0.142107 -0.142106 0.979598	0.274278 -0.274277 0.921707	0.322516 -0.215498 0.921707	0.142107 -0.142106 0.979598	0.322516 -0.215498 0.921707	0.1671 -0.111652 0.979598	0.1671 -0.111652 0.979598	0.322516 -0.215498 0.921707	0.358361 -0.148437 0.921707	0.1671 -0.111652 0.979598	0.358361 -0.148437 0.921707	0.185671 -0.0769072 0.979598	0.185671 -0.0769072 0.979598	0.358361 -0.148437 0.921707	0.380434 -0.0756724 0.921707	0.185671 -0.0769072 0.979598	0.380434 -0.0756724 0.921707	0.197107 -0.0392068 0.979598	0.197107 -0.0392068 0.979598	0.380434 -0.0756724 0.921707	0.387887 5.70723e-007 0.921707	0.197107 -0.0392068 0.979598	0.387887 5.70723e-007 0.921707	0.200969 3.22768e-007 0.979598	0.200969 3.22768e-007 0.979598	0.387887 5.70723e-007 0.921707	0.380434 0.0756737 0.921707	0.200969 3.22768e-007 0.979598	0.380434 0.0756737 0.921707	0.197107 0.0392074 0.979598	0.197107 0.0392074 0.979598	0.380434 0.0756737 0.921707	0.358361 0.148439 0.921707	0.197107 0.0392074 0.979598	0.358361 0.148439 0.921707	0.185671 0.0769077 0.979598	0.185671 0.0769077 0.979598	0.358361 0.148439 0.921707	0.322516 0.215499 0.921707	0.185671 0.0769077 0.979598	0.322516 0.215499 0.921707	0.167099 0.111653 0.979598	0.167099 0.111653 0.979598	0.322516 0.215499 0.921707	0.274277 0.274278 0.921707	0.167099 0.111653 0.979598	0.274277 0.274278 0.921707	0.142106 0.142107 0.979598	0.142106 0.142107 0.979598	0.274277 0.274278 0.921707	0.215498 0.322517 0.921707	0.142106 0.142107 0.979598	0.215498 0.322517 0.921707	0.111652 0.1671 0.979598	0.111652 0.1671 0.979598	0.215498 0.322517 0.921707	0.148437 0.358361 0.921707	0.111652 0.1671 0.979598	0.148437 0.358361 0.921707	0.0769071 0.185671 0.979598	0.0769071 0.185671 0.979598	0.148437 0.358361 0.921707	0.0756724 0.380434 0.921707	0.0769071 0.185671 0.979598	0.0756724 0.380434 0.921707	0.0392067 0.197107 0.979598	0.0392067 0.197107 0.979598	0.0756724 0.380434 0.921707	-2.84465e-007 0.387887 0.921707	0.0392067 0.197107 0.979598	-2.84465e-007 0.387887 0.921707	0.0 0.200969 0.979598	-2.84465e-007 0.387887 0.921707	-3.50847e-007 0.559775 0.828645	-0.109207 0.549019 0.828645	-2.84465e-007 0.387887 0.921707	-0.109207 0.549019 0.828645	-0.075673 0.380434 0.921707	-0.075673 0.380434 0.921707	-0.109207 0.549019 0.828645	-0.214217 0.517165 0.828645	-0.075673 0.380434 0.921707	-0.214217 0.517165 0.828645	-0.148438 0.358361 0.921707	-0.148438 0.358361 0.921707	-0.214217 0.517165 0.828645	-0.310994 0.465436 0.828645	-0.148438 0.358361 0.921707	-0.310994 0.465436 0.828645	-0.215498 0.322516 0.921707	-0.215498 0.322516 0.921707	-0.310994 0.465436 0.828645	-0.395821 0.395821 0.828645	-0.215498 0.322516 0.921707	-0.395821 0.395821 0.828645	-0.274277 0.274277 0.921707	-0.274277 0.274277 0.921707	-0.395821 0.395821 0.828645	-0.465436 0.310994 0.828645	-0.274277 0.274277 0.921707	-0.465436 0.310994 0.828645	-0.322516 0.215498 0.921707	-0.322516 0.215498 0.921707	-0.465436 0.310994 0.828645	-0.517165 0.214217 0.828645	-0.322516 0.215498 0.921707	-0.517165 0.214217 0.828645	-0.358361 0.148438 0.921707	-0.358361 0.148438 0.921707	-0.517165 0.214217 0.828645	-0.549019 0.109207 0.828645	-0.358361 0.148438 0.921707	-0.549019 0.109207 0.828645	-0.380434 0.0756729 0.921707	-0.380434 0.0756729 0.921707	-0.549019 0.109207 0.828645	-0.559775 -1.60456e-007 0.828645	-0.380434 0.0756729 0.921707	-0.559775 -1.60456e-007 0.828645	-0.387887 0.0 0.921707	-0.387887 0.0 0.921707	-0.559775 -1.60456e-007 0.828645	-0.549019 -0.109207 0.828645	-0.387887 0.0 0.921707	-0.549019 -0.109207 0.828645	-0.380434 -0.0756732 0.921707	-0.380434 -0.0756732 0.921707	-0.549019 -0.109207 0.828645	-0.517165 -0.214217 0.828645	-0.380434 -0.0756732 0.921707	-0.517165 -0.214217 0.828645	-0.358361 -0.148438 0.921707	-0.358361 -0.148438 0.921707	-0.517165 -0.214217 0.828645	-0.465436 -0.310995 0.828645	-0.358361 -0.148438 0.921707	-0.465436 -0.310995 0.828645	-0.322516 -0.215499 0.921707	-0.322516 -0.215499 0.921707	-0.465436 -0.310995 0.828645	-0.39582 -0.395821 0.828645	-0.322516 -0.215499 0.921707	-0.39582 -0.395821 0.828645	-0.274277 -0.274278 0.921707	-0.274277 -0.274278 0.921707	-0.39582 -0.395821 0.828645	-0.310994 -0.465436 0.828645	-0.274277 -0.274278 0.921707	-0.310994 -0.465436 0.828645	-0.215498 -0.322516 0.921707	-0.215498 -0.322516 0.921707	-0.310994 -0.465436 0.828645	-0.214216 -0.517165 0.828645	-0.215498 -0.322516 0.921707	-0.214216 -0.517165 0.828645	-0.148438 -0.358361 0.921707	-0.148438 -0.358361 0.921707	-0.214216 -0.517165 0.828645	-0.109206 -0.549019 0.828645	-0.148438 -0.358361 0.921707	-0.109206 -0.549019 0.828645	-0.0756727 -0.380434 0.921707	-0.0756727 -0.380434 0.921707	-0.109206 -0.549019 0.828645	4.65801e-007 -0.559775 0.828645	-0.0756727 -0.380434 0.921707	4.65801e-007 -0.559775 0.828645	3.96817e-007 -0.387887 0.921707	3.96817e-007 -0.387887 0.921707	4.65801e-007 -0.559775 0.828645	0.109207 -0.549019 0.828645	3.96817e-007 -0.387887 0.921707	0.109207 -0.549019 0.828645	0.0756734 -0.380434 0.921707	0.0756734 -0.380434 0.921707	0.109207 -0.549019 0.828645	0.214217 -0.517164 0.828645	0.0756734 -0.380434 0.921707	0.214217 -0.517164 0.828645	0.148438 -0.358361 0.921707	0.148438 -0.358361 0.921707	0.214217 -0.517164 0.828645	0.310995 -0.465436 0.828645	0.148438 -0.358361 0.921707	0.310995 -0.465436 0.828645	0.215499 -0.322516 0.921707	0.215499 -0.322516 0.921707	0.310995 -0.465436 0.828645	0.395821 -0.39582 0.828645	0.215499 -0.322516 0.921707	0.395821 -0.39582 0.828645	0.274278 -0.274277 0.921707	0.274278 -0.274277 0.921707	0.395821 -0.39582 0.828645	0.465436 -0.310994 0.828645	0.274278 -0.274277 0.921707	0.465436 -0.310994 0.828645	0.322516 -0.215498 0.921707	0.322516 -0.215498 0.921707	0.465436 -0.310994 0.828645	0.517165 -0.214216 0.828645	0.322516 -0.215498 0.921707	0.517165 -0.214216 0.828645	0.358361 -0.148437 0.921707	0.358361 -0.148437 0.921707	0.517165 -0.214216 0.828645	0.549019 -0.109206 0.828645	0.358361 -0.148437 0.921707	0.549019 -0.109206 0.828645	0.380434 -0.0756724 0.921707	0.380434 -0.0756724 0.921707	0.549019 -0.109206 0.828645	0.559775 7.927e-007 0.828645	0.380434 -0.0756724 0.921707	0.559775 7.927e-007 0.828645	0.387887 5.70723e-007 0.921707	0.387887 5.70723e-007 0.921707	0.559775 7.927e-007 0.828645	0.549019 0.109208 0.828645	0.387887 5.70723e-007 0.921707	0.549019 0.109208 0.828645	0.380434 0.0756737 0.921707	0.380434 0.0756737 0.921707	0.549019 0.109208 0.828645	0.517164 0.214218 0.828645	0.380434 0.0756737 0.921707	0.517164 0.214218 0.828645	0.358361 0.148439 0.921707	0.358361 0.148439 0.921707	0.517164 0.214218 0.828645	0.465435 0.310995 0.828645	0.358361 0.148439 0.921707	0.465435 0.310995 0.828645	0.322516 0.215499 0.921707	0.322516 0.215499 0.921707	0.465435 0.310995 0.828645	0.39582 0.395821 0.828645	0.322516 0.215499 0.921707	0.39582 0.395821 0.828645	0.274277 0.274278 0.921707	0.274277 0.274278 0.921707	0.39582 0.395821 0.828645	0.310993 0.465437 0.828645	0.274277 0.274278 0.921707	0.310993 0.465437 0.828645	0.215498 0.322517 0.921707	0.215498 0.322517 0.921707	0.310993 0.465437 0.828645	0.214216 0.517165 0.828645	0.215498 0.322517 0.921707	0.214216 0.517165 0.828645	0.148437 0.358361 0.921707	0.148437 0.358361 0.921707	0.214216 0.517165 0.828645	0.109206 0.549019 0.828645	0.148437 0.358361 0.921707	0.109206 0.549019 0.828645	0.0756724 0.380434 0.921707	0.0756724 0.380434 0.921707	0.109206 0.549019 0.828645	-3.50847e-007 0.559775 0.828645	0.0756724 0.380434 0.921707	-3.50847e-007 0.559775 0.828645	-2.84465e-007 0.387887 0.921707	-3.50847e-007 0.559775 0.828645	-4.37867e-007 0.710141 0.704059	-0.138542 0.696496 0.704059	-3.50847e-007 0.559775 0.828645	-0.138542 0.696496 0.704059	-0.109207 0.549019 0.828645	-0.109207 0.549019 0.828645	-0.138542 0.696496 0.704059	-0.271759 0.656085 0.704059	-0.109207 0.549019 0.828645	-0.271759 0.656085 0.704059	-0.214217 0.517165 0.828645	-0.214217 0.517165 0.828645	-0.271759 0.656085 0.704059	-0.394533 0.590461 0.704059	-0.214217 0.517165 0.828645	-0.394533 0.590461 0.704059	-0.310994 0.465436 0.828645	-0.310994 0.465436 0.828645	-0.394533 0.590461 0.704059	-0.502146 0.502146 0.704059	-0.310994 0.465436 0.828645	-0.502146 0.502146 0.704059	-0.395821 0.395821 0.828645	-0.395821 0.395821 0.828645	-0.502146 0.502146 0.704059	-0.590461 0.394533 0.704059	-0.395821 0.395821 0.828645	-0.590461 0.394533 0.704059	-0.465436 0.310994 0.828645	-0.465436 0.310994 0.828645	-0.590461 0.394533 0.704059	-0.656085 0.271759 0.704059	-0.465436 0.310994 0.828645	-0.656085 0.271759 0.704059	-0.517165 0.214217 0.828645	-0.517165 0.214217 0.828645	-0.656085 0.271759 0.704059	-0.696496 0.138541 0.704059	-0.517165 0.214217 0.828645	-0.696496 0.138541 0.704059	-0.549019 0.109207 0.828645	-0.549019 0.109207 0.828645	-0.696496 0.138541 0.704059	-0.710141 -2.9511e-007 0.704059	-0.549019 0.109207 0.828645	-0.710141 -2.9511e-007 0.704059	-0.559775 -1.60456e-007 0.828645	-0.559775 -1.60456e-007 0.828645	-0.710141 -2.9511e-007 0.704059	-0.696496 -0.138542 0.704059	-0.559775 -1.60456e-007 0.828645	-0.696496 -0.138542 0.704059	-0.549019 -0.109207 0.828645	-0.549019 -0.109207 0.828645	-0.696496 -0.138542 0.704059	-0.656085 -0.27176 0.704059	-0.549019 -0.109207 0.828645	-0.656085 -0.27176 0.704059	-0.517165 -0.214217 0.828645	-0.517165 -0.214217 0.828645	-0.656085 -0.27176 0.704059	-0.590461 -0.394534 0.704059	-0.517165 -0.214217 0.828645	-0.590461 -0.394534 0.704059	-0.465436 -0.310995 0.828645	-0.465436 -0.310995 0.828645	-0.590461 -0.394534 0.704059	-0.502145 -0.502146 0.704059	-0.465436 -0.310995 0.828645	-0.502145 -0.502146 0.704059	-0.39582 -0.395821 0.828645	-0.39582 -0.395821 0.828645	-0.502145 -0.502146 0.704059	-0.394533 -0.590461 0.704059	-0.39582 -0.395821 0.828645	-0.394533 -0.590461 0.704059	-0.310994 -0.465436 0.828645	-0.310994 -0.465436 0.828645	-0.394533 -0.590461 0.704059	-0.271759 -0.656085 0.704059	-0.310994 -0.465436 0.828645	-0.271759 -0.656085 0.704059	-0.214216 -0.517165 0.828645	-0.214216 -0.517165 0.828645	-0.271759 -0.656085 0.704059	-0.138541 -0.696496 0.704059	-0.214216 -0.517165 0.828645	-0.138541 -0.696496 0.704059	-0.109206 -0.549019 0.828645	-0.109206 -0.549019 0.828645	-0.138541 -0.696496 0.704059	6.17812e-007 -0.710141 0.704059	-0.109206 -0.549019 0.828645	6.17812e-007 -0.710141 0.704059	4.65801e-007 -0.559775 0.828645	4.65801e-007 -0.559775 0.828645	6.17812e-007 -0.710141 0.704059	0.138542 -0.696496 0.704059	4.65801e-007 -0.559775 0.828645	0.138542 -0.696496 0.704059	0.109207 -0.549019 0.828645	0.109207 -0.549019 0.828645	0.138542 -0.696496 0.704059	0.27176 -0.656085 0.704059	0.109207 -0.549019 0.828645	0.27176 -0.656085 0.704059	0.214217 -0.517164 0.828645	0.214217 -0.517164 0.828645	0.27176 -0.656085 0.704059	0.394534 -0.590461 0.704059	0.214217 -0.517164 0.828645	0.394534 -0.590461 0.704059	0.310995 -0.465436 0.828645	0.310995 -0.465436 0.828645	0.394534 -0.590461 0.704059	0.502146 -0.502145 0.704059	0.310995 -0.465436 0.828645	0.502146 -0.502145 0.704059	0.395821 -0.39582 0.828645	0.395821 -0.39582 0.828645	0.502146 -0.502145 0.704059	0.590461 -0.394533 0.704059	0.395821 -0.39582 0.828645	0.590461 -0.394533 0.704059	0.465436 -0.310994 0.828645	0.465436 -0.310994 0.828645	0.590461 -0.394533 0.704059	0.656086 -0.271758 0.704059	0.465436 -0.310994 0.828645	0.656086 -0.271758 0.704059	0.517165 -0.214216 0.828645	0.517165 -0.214216 0.828645	0.656086 -0.271758 0.704059	0.696496 -0.138541 0.704059	0.517165 -0.214216 0.828645	0.696496 -0.138541 0.704059	0.549019 -0.109206 0.828645	0.549019 -0.109206 0.828645	0.696496 -0.138541 0.704059	0.710141 1.09647e-006 0.704059	0.549019 -0.109206 0.828645	0.710141 1.09647e-006 0.704059	0.559775 7.927e-007 0.828645	0.559775 7.927e-007 0.828645	0.710141 1.09647e-006 0.704059	0.696496 0.138543 0.704059	0.559775 7.927e-007 0.828645	0.696496 0.138543 0.704059	0.549019 0.109208 0.828645	0.549019 0.109208 0.828645	0.696496 0.138543 0.704059	0.656085 0.27176 0.704059	0.549019 0.109208 0.828645	0.656085 0.27176 0.704059	0.517164 0.214218 0.828645	0.517164 0.214218 0.828645	0.656085 0.27176 0.704059	0.59046 0.394535 0.704059	0.517164 0.214218 0.828645	0.59046 0.394535 0.704059	0.465435 0.310995 0.828645	0.465435 0.310995 0.828645	0.59046 0.394535 0.704059	0.502145 0.502147 0.704059	0.465435 0.310995 0.828645	0.502145 0.502147 0.704059	0.39582 0.395821 0.828645	0.39582 0.395821 0.828645	0.502145 0.502147 0.704059	0.394532 0.590462 0.704059	0.39582 0.395821 0.828645	0.394532 0.590462 0.704059	0.310993 0.465437 0.828645	0.310993 0.465437 0.828645	0.394532 0.590462 0.704059	0.271758 0.656086 0.704059	0.310993 0.465437 0.828645	0.271758 0.656086 0.704059	0.214216 0.517165 0.828645	0.214216 0.517165 0.828645	0.271758 0.656086 0.704059	0.138541 0.696496 0.704059	0.214216 0.517165 0.828645	0.138541 0.696496 0.704059	0.109206 0.549019 0.828645	0.109206 0.549019 0.828645	0.138541 0.696496 0.704059	-4.37867e-007 0.710141 0.704059	0.109206 0.549019 0.828645	-4.37867e-007 0.710141 0.704059	-3.50847e-007 0.559775 0.828645	-4.37867e-007 0.710141 0.704059	-5.07089e-007 0.83334 0.552761	-0.162577 0.817327 0.552761	-4.37867e-007 0.710141 0.704059	-0.162577 0.817327 0.552761	-0.138542 0.696496 0.704059	-0.138542 0.696496 0.704059	-0.162577 0.817327 0.552761	-0.318905 0.769906 0.552761	-0.138542 0.696496 0.704059	-0.318905 0.769906 0.552761	-0.271759 0.656085 0.704059	-0.271759 0.656085 0.704059	-0.318905 0.769906 0.552761	-0.462979 0.692897 0.552761	-0.271759 0.656085 0.704059	-0.462979 0.692897 0.552761	-0.394533 0.590461 0.704059	-0.394533 0.590461 0.704059	-0.462979 0.692897 0.552761	-0.58926 0.58926 0.552761	-0.394533 0.590461 0.704059	-0.58926 0.58926 0.552761	-0.502146 0.502146 0.704059	-0.502146 0.502146 0.704059	-0.58926 0.58926 0.552761	-0.692897 0.462979 0.552761	-0.502146 0.502146 0.704059	-0.692897 0.462979 0.552761	-0.590461 0.394533 0.704059	-0.590461 0.394533 0.704059	-0.692897 0.462979 0.552761	-0.769906 0.318905 0.552761	-0.590461 0.394533 0.704059	-0.769906 0.318905 0.552761	-0.656085 0.271759 0.704059	-0.656085 0.271759 0.704059	-0.769906 0.318905 0.552761	-0.817327 0.162576 0.552761	-0.656085 0.271759 0.704059	-0.817327 0.162576 0.552761	-0.696496 0.138541 0.704059	-0.696496 0.138541 0.704059	-0.817327 0.162576 0.552761	-0.83334 -3.0882e-007 0.552761	-0.696496 0.138541 0.704059	-0.83334 -3.0882e-007 0.552761	-0.710141 -2.9511e-007 0.704059	-0.710141 -2.9511e-007 0.704059	-0.83334 -3.0882e-007 0.552761	-0.817327 -0.162577 0.552761	-0.710141 -2.9511e-007 0.704059	-0.817327 -0.162577 0.552761	-0.696496 -0.138542 0.704059	-0.696496 -0.138542 0.704059	-0.817327 -0.162577 0.552761	-0.769905 -0.318906 0.552761	-0.696496 -0.138542 0.704059	-0.769905 -0.318906 0.552761	-0.656085 -0.27176 0.704059	-0.656085 -0.27176 0.704059	-0.769905 -0.318906 0.552761	-0.692896 -0.462979 0.552761	-0.656085 -0.27176 0.704059	-0.692896 -0.462979 0.552761	-0.590461 -0.394534 0.704059	-0.590461 -0.394534 0.704059	-0.692896 -0.462979 0.552761	-0.58926 -0.589261 0.552761	-0.590461 -0.394534 0.704059	-0.58926 -0.589261 0.552761	-0.502145 -0.502146 0.704059	-0.502145 -0.502146 0.704059	-0.58926 -0.589261 0.552761	-0.462978 -0.692897 0.552761	-0.502145 -0.502146 0.704059	-0.462978 -0.692897 0.552761	-0.394533 -0.590461 0.704059	-0.394533 -0.590461 0.704059	-0.462978 -0.692897 0.552761	-0.318905 -0.769906 0.552761	-0.394533 -0.590461 0.704059	-0.318905 -0.769906 0.552761	-0.271759 -0.656085 0.704059	-0.271759 -0.656085 0.704059	-0.318905 -0.769906 0.552761	-0.162576 -0.817327 0.552761	-0.271759 -0.656085 0.704059	-0.162576 -0.817327 0.552761	-0.138541 -0.696496 0.704059	-0.138541 -0.696496 0.704059	-0.162576 -0.817327 0.552761	6.69309e-007 -0.83334 0.552761	-0.138541 -0.696496 0.704059	6.69309e-007 -0.83334 0.552761	6.17812e-007 -0.710141 0.704059	6.17812e-007 -0.710141 0.704059	6.69309e-007 -0.83334 0.552761	0.162577 -0.817327 0.552761	6.17812e-007 -0.710141 0.704059	0.162577 -0.817327 0.552761	0.138542 -0.696496 0.704059	0.138542 -0.696496 0.704059	0.162577 -0.817327 0.552761	0.318906 -0.769905 0.552761	0.138542 -0.696496 0.704059	0.318906 -0.769905 0.552761	0.27176 -0.656085 0.704059	0.27176 -0.656085 0.704059	0.318906 -0.769905 0.552761	0.46298 -0.692896 0.552761	0.27176 -0.656085 0.704059	0.46298 -0.692896 0.552761	0.394534 -0.590461 0.704059	0.394534 -0.590461 0.704059	0.46298 -0.692896 0.552761	0.589261 -0.589259 0.552761	0.394534 -0.590461 0.704059	0.589261 -0.589259 0.552761	0.502146 -0.502145 0.704059	0.502146 -0.502145 0.704059	0.589261 -0.589259 0.552761	0.692897 -0.462978 0.552761	0.502146 -0.502145 0.704059	0.692897 -0.462978 0.552761	0.590461 -0.394533 0.704059	0.590461 -0.394533 0.704059	0.692897 -0.462978 0.552761	0.769906 -0.318904 0.552761	0.590461 -0.394533 0.704059	0.769906 -0.318904 0.552761	0.656086 -0.271758 0.704059	0.656086 -0.271758 0.704059	0.769906 -0.318904 0.552761	0.817328 -0.162575 0.552761	0.656086 -0.271758 0.704059	0.817328 -0.162575 0.552761	0.696496 -0.138541 0.704059	0.696496 -0.138541 0.704059	0.817328 -0.162575 0.552761	0.83334 1.27974e-006 0.552761	0.696496 -0.138541 0.704059	0.83334 1.27974e-006 0.552761	0.710141 1.09647e-006 0.704059	0.710141 1.09647e-006 0.704059	0.83334 1.27974e-006 0.552761	0.817327 0.162578 0.552761	0.710141 1.09647e-006 0.704059	0.817327 0.162578 0.552761	0.696496 0.138543 0.704059	0.696496 0.138543 0.704059	0.817327 0.162578 0.552761	0.769905 0.318907 0.552761	0.696496 0.138543 0.704059	0.769905 0.318907 0.552761	0.656085 0.27176 0.704059	0.656085 0.27176 0.704059	0.769905 0.318907 0.552761	0.692896 0.46298 0.552761	0.656085 0.27176 0.704059	0.692896 0.46298 0.552761	0.59046 0.394535 0.704059	0.59046 0.394535 0.704059	0.692896 0.46298 0.552761	0.589259 0.589261 0.552761	0.59046 0.394535 0.704059	0.589259 0.589261 0.552761	0.502145 0.502147 0.704059	0.502145 0.502147 0.704059	0.589259 0.589261 0.552761	0.462977 0.692897 0.552761	0.502145 0.502147 0.704059	0.462977 0.692897 0.552761	0.394532 0.590462 0.704059	0.394532 0.590462 0.704059	0.462977 0.692897 0.552761	0.318904 0.769906 0.552761	0.394532 0.590462 0.704059	0.318904 0.769906 0.552761	0.271758 0.656086 0.704059	0.271758 0.656086 0.704059	0.318904 0.769906 0.552761	0.162575 0.817328 0.552761	0.271758 0.656086 0.704059	0.162575 0.817328 0.552761	0.138541 0.696496 0.704059	0.138541 0.696496 0.704059	0.162575 0.817328 0.552761	-5.07089e-007 0.83334 0.552761	0.138541 0.696496 0.704059	-5.07089e-007 0.83334 0.552761	-4.37867e-007 0.710141 0.704059	-5.07089e-007 0.83334 0.552761	-6.35304e-007 0.924766 0.380537	-0.180413 0.906997 0.380537	-5.07089e-007 0.83334 0.552761	-0.180413 0.906997 0.380537	-0.162577 0.817327 0.552761	-0.162577 0.817327 0.552761	-0.180413 0.906997 0.380537	-0.353892 0.854372 0.380537	-0.162577 0.817327 0.552761	-0.353892 0.854372 0.380537	-0.318905 0.769906 0.552761	-0.318905 0.769906 0.552761	-0.353892 0.854372 0.380537	-0.513772 0.768915 0.380537	-0.318905 0.769906 0.552761	-0.513772 0.768915 0.380537	-0.462979 0.692897 0.552761	-0.462979 0.692897 0.552761	-0.513772 0.768915 0.380537	-0.653908 0.653908 0.380537	-0.462979 0.692897 0.552761	-0.653908 0.653908 0.380537	-0.58926 0.58926 0.552761	-0.58926 0.58926 0.552761	-0.653908 0.653908 0.380537	-0.768915 0.513772 0.380537	-0.58926 0.58926 0.552761	-0.768915 0.513772 0.380537	-0.692897 0.462979 0.552761	-0.692897 0.462979 0.552761	-0.768915 0.513772 0.380537	-0.854372 0.353892 0.380537	-0.692897 0.462979 0.552761	-0.854372 0.353892 0.380537	-0.769906 0.318905 0.552761	-0.769906 0.318905 0.552761	-0.854372 0.353892 0.380537	-0.906997 0.180413 0.380537	-0.769906 0.318905 0.552761	-0.906997 0.180413 0.380537	-0.817327 0.162576 0.552761	-0.817327 0.162576 0.552761	-0.906997 0.180413 0.380537	-0.924766 -2.03345e-007 0.380537	-0.817327 0.162576 0.552761	-0.924766 -2.03345e-007 0.380537	-0.83334 -3.0882e-007 0.552761	-0.83334 -3.0882e-007 0.552761	-0.924766 -2.03345e-007 0.380537	-0.906996 -0.180413 0.380537	-0.83334 -3.0882e-007 0.552761	-0.906996 -0.180413 0.380537	-0.817327 -0.162577 0.552761	-0.817327 -0.162577 0.552761	-0.906996 -0.180413 0.380537	-0.854372 -0.353893 0.380537	-0.817327 -0.162577 0.552761	-0.854372 -0.353893 0.380537	-0.769905 -0.318906 0.552761	-0.769905 -0.318906 0.552761	-0.854372 -0.353893 0.380537	-0.768914 -0.513773 0.380537	-0.769905 -0.318906 0.552761	-0.768914 -0.513773 0.380537	-0.692896 -0.462979 0.552761	-0.692896 -0.462979 0.552761	-0.768914 -0.513773 0.380537	-0.653908 -0.653908 0.380537	-0.692896 -0.462979 0.552761	-0.653908 -0.653908 0.380537	-0.58926 -0.589261 0.552761	-0.58926 -0.589261 0.552761	-0.653908 -0.653908 0.380537	-0.513772 -0.768915 0.380537	-0.58926 -0.589261 0.552761	-0.513772 -0.768915 0.380537	-0.462978 -0.692897 0.552761	-0.462978 -0.692897 0.552761	-0.513772 -0.768915 0.380537	-0.353892 -0.854372 0.380537	-0.462978 -0.692897 0.552761	-0.353892 -0.854372 0.380537	-0.318905 -0.769906 0.552761	-0.318905 -0.769906 0.552761	-0.353892 -0.854372 0.380537	-0.180412 -0.906997 0.380537	-0.318905 -0.769906 0.552761	-0.180412 -0.906997 0.380537	-0.162576 -0.817327 0.552761	-0.162576 -0.817327 0.552761	-0.180412 -0.906997 0.380537	8.07366e-007 -0.924766 0.380537	-0.162576 -0.817327 0.552761	8.07366e-007 -0.924766 0.380537	6.69309e-007 -0.83334 0.552761	6.69309e-007 -0.83334 0.552761	8.07366e-007 -0.924766 0.380537	0.180414 -0.906996 0.380537	6.69309e-007 -0.83334 0.552761	0.180414 -0.906996 0.380537	0.162577 -0.817327 0.552761	0.162577 -0.817327 0.552761	0.180414 -0.906996 0.380537	0.353893 -0.854372 0.380537	0.162577 -0.817327 0.552761	0.353893 -0.854372 0.380537	0.318906 -0.769905 0.552761	0.318906 -0.769905 0.552761	0.353893 -0.854372 0.380537	0.513773 -0.768914 0.380537	0.318906 -0.769905 0.552761	0.513773 -0.768914 0.380537	0.46298 -0.692896 0.552761	0.46298 -0.692896 0.552761	0.513773 -0.768914 0.380537	0.653909 -0.653907 0.380537	0.46298 -0.692896 0.552761	0.653909 -0.653907 0.380537	0.589261 -0.589259 0.552761	0.589261 -0.589259 0.552761	0.653909 -0.653907 0.380537	0.768915 -0.513771 0.380537	0.589261 -0.589259 0.552761	0.768915 -0.513771 0.380537	0.692897 -0.462978 0.552761	0.692897 -0.462978 0.552761	0.768915 -0.513771 0.380537	0.854373 -0.353891 0.380537	0.692897 -0.462978 0.552761	0.854373 -0.353891 0.380537	0.769906 -0.318904 0.552761	0.769906 -0.318904 0.552761	0.854373 -0.353891 0.380537	0.906997 -0.180411 0.380537	0.769906 -0.318904 0.552761	0.906997 -0.180411 0.380537	0.817328 -0.162575 0.552761	0.817328 -0.162575 0.552761	0.906997 -0.180411 0.380537	0.924766 1.44508e-006 0.380537	0.817328 -0.162575 0.552761	0.924766 1.44508e-006 0.380537	0.83334 1.27974e-006 0.552761	0.83334 1.27974e-006 0.552761	0.924766 1.44508e-006 0.380537	0.906996 0.180414 0.380537	0.83334 1.27974e-006 0.552761	0.906996 0.180414 0.380537	0.817327 0.162578 0.552761	0.817327 0.162578 0.552761	0.906996 0.180414 0.380537	0.854371 0.353894 0.380537	0.817327 0.162578 0.552761	0.854371 0.353894 0.380537	0.769905 0.318907 0.552761	0.769905 0.318907 0.552761	0.854371 0.353894 0.380537	0.768914 0.513774 0.380537	0.769905 0.318907 0.552761	0.768914 0.513774 0.380537	0.692896 0.46298 0.552761	0.692896 0.46298 0.552761	0.768914 0.513774 0.380537	0.653907 0.653909 0.380537	0.692896 0.46298 0.552761	0.653907 0.653909 0.380537	0.589259 0.589261 0.552761	0.589259 0.589261 0.552761	0.653907 0.653909 0.380537	0.513771 0.768915 0.380537	0.589259 0.589261 0.552761	0.513771 0.768915 0.380537	0.462977 0.692897 0.552761	0.462977 0.692897 0.552761	0.513771 0.768915 0.380537	0.353891 0.854373 0.380537	0.462977 0.692897 0.552761	0.353891 0.854373 0.380537	0.318904 0.769906 0.552761	0.318904 0.769906 0.552761	0.353891 0.854373 0.380537	0.180411 0.906997 0.380537	0.318904 0.769906 0.552761	0.180411 0.906997 0.380537	0.162575 0.817328 0.552761	0.162575 0.817328 0.552761	0.180411 0.906997 0.380537	-6.35304e-007 0.924766 0.380537	0.162575 0.817328 0.552761	-6.35304e-007 0.924766 0.380537	-5.07089e-007 0.83334 0.552761	-6.35304e-007 0.924766 0.380537	-7.6109e-007 0.981015 0.19393	-0.191387 0.962165 0.19393	-6.35304e-007 0.924766 0.380537	-0.191387 0.962165 0.19393	-0.180413 0.906997 0.380537	-0.180413 0.906997 0.380537	-0.191387 0.962165 0.19393	-0.375418 0.90634 0.19393	-0.180413 0.906997 0.380537	-0.375418 0.90634 0.19393	-0.353892 0.854372 0.380537	-0.353892 0.854372 0.380537	-0.375418 0.90634 0.19393	-0.545023 0.815684 0.19393	-0.353892 0.854372 0.380537	-0.545023 0.815684 0.19393	-0.513772 0.768915 0.380537	-0.513772 0.768915 0.380537	-0.545023 0.815684 0.19393	-0.693683 0.693683 0.19393	-0.513772 0.768915 0.380537	-0.693683 0.693683 0.19393	-0.653908 0.653908 0.380537	-0.653908 0.653908 0.380537	-0.693683 0.693683 0.19393	-0.815684 0.545023 0.19393	-0.653908 0.653908 0.380537	-0.815684 0.545023 0.19393	-0.768915 0.513772 0.380537	-0.768915 0.513772 0.380537	-0.815684 0.545023 0.19393	-0.90634 0.375418 0.19393	-0.768915 0.513772 0.380537	-0.90634 0.375418 0.19393	-0.854372 0.353892 0.380537	-0.854372 0.353892 0.380537	-0.90634 0.375418 0.19393	-0.962165 0.191386 0.19393	-0.854372 0.353892 0.380537	-0.962165 0.191386 0.19393	-0.906997 0.180413 0.380537	-0.906997 0.180413 0.380537	-0.962165 0.191386 0.19393	-0.981015 -3.15515e-007 0.19393	-0.906997 0.180413 0.380537	-0.981015 -3.15515e-007 0.19393	-0.924766 -2.03345e-007 0.380537	-0.924766 -2.03345e-007 0.380537	-0.981015 -3.15515e-007 0.19393	-0.962165 -0.191387 0.19393	-0.924766 -2.03345e-007 0.380537	-0.962165 -0.191387 0.19393	-0.906996 -0.180413 0.380537	-0.906996 -0.180413 0.380537	-0.962165 -0.191387 0.19393	-0.90634 -0.375419 0.19393	-0.906996 -0.180413 0.380537	-0.90634 -0.375419 0.19393	-0.854372 -0.353893 0.380537	-0.854372 -0.353893 0.380537	-0.90634 -0.375419 0.19393	-0.815684 -0.545023 0.19393	-0.854372 -0.353893 0.380537	-0.815684 -0.545023 0.19393	-0.768914 -0.513773 0.380537	-0.768914 -0.513773 0.380537	-0.815684 -0.545023 0.19393	-0.693682 -0.693683 0.19393	-0.768914 -0.513773 0.380537	-0.693682 -0.693683 0.19393	-0.653908 -0.653908 0.380537	-0.653908 -0.653908 0.380537	-0.693682 -0.693683 0.19393	-0.545022 -0.815685 0.19393	-0.653908 -0.653908 0.380537	-0.545022 -0.815685 0.19393	-0.513772 -0.768915 0.380537	-0.513772 -0.768915 0.380537	-0.545022 -0.815685 0.19393	-0.375418 -0.90634 0.19393	-0.513772 -0.768915 0.380537	-0.375418 -0.90634 0.19393	-0.353892 -0.854372 0.380537	-0.353892 -0.854372 0.380537	-0.375418 -0.90634 0.19393	-0.191386 -0.962166 0.19393	-0.353892 -0.854372 0.380537	-0.191386 -0.962166 0.19393	-0.180412 -0.906997 0.380537	-0.180412 -0.906997 0.380537	-0.191386 -0.962166 0.19393	9.51363e-007 -0.981015 0.19393	-0.180412 -0.906997 0.380537	9.51363e-007 -0.981015 0.19393	8.07366e-007 -0.924766 0.380537	8.07366e-007 -0.924766 0.380537	9.51363e-007 -0.981015 0.19393	0.191388 -0.962165 0.19393	8.07366e-007 -0.924766 0.380537	0.191388 -0.962165 0.19393	0.180414 -0.906996 0.380537	0.180414 -0.906996 0.380537	0.191388 -0.962165 0.19393	0.375419 -0.90634 0.19393	0.180414 -0.906996 0.380537	0.375419 -0.90634 0.19393	0.353893 -0.854372 0.380537	0.353893 -0.854372 0.380537	0.375419 -0.90634 0.19393	0.545024 -0.815684 0.19393	0.353893 -0.854372 0.380537	0.545024 -0.815684 0.19393	0.513773 -0.768914 0.380537	0.513773 -0.768914 0.380537	0.545024 -0.815684 0.19393	0.693684 -0.693682 0.19393	0.513773 -0.768914 0.380537	0.693684 -0.693682 0.19393	0.653909 -0.653907 0.380537	0.653909 -0.653907 0.380537	0.693684 -0.693682 0.19393	0.815685 -0.545022 0.19393	0.653909 -0.653907 0.380537	0.815685 -0.545022 0.19393	0.768915 -0.513771 0.380537	0.768915 -0.513771 0.380537	0.815685 -0.545022 0.19393	0.906341 -0.375417 0.19393	0.768915 -0.513771 0.380537	0.906341 -0.375417 0.19393	0.854373 -0.353891 0.380537	0.854373 -0.353891 0.380537	0.906341 -0.375417 0.19393	0.962166 -0.191385 0.19393	0.854373 -0.353891 0.380537	0.962166 -0.191385 0.19393	0.906997 -0.180411 0.380537	0.906997 -0.180411 0.380537	0.962166 -0.191385 0.19393	0.981015 1.57758e-006 0.19393	0.906997 -0.180411 0.380537	0.981015 1.57758e-006 0.19393	0.924766 1.44508e-006 0.380537	0.924766 1.44508e-006 0.380537	0.981015 1.57758e-006 0.19393	0.962165 0.191388 0.19393	0.924766 1.44508e-006 0.380537	0.962165 0.191388 0.19393	0.906996 0.180414 0.380537	0.906996 0.180414 0.380537	0.962165 0.191388 0.19393	0.906339 0.37542 0.19393	0.906996 0.180414 0.380537	0.906339 0.37542 0.19393	0.854371 0.353894 0.380537	0.854371 0.353894 0.380537	0.906339 0.37542 0.19393	0.815683 0.545024 0.19393	0.854371 0.353894 0.380537	0.815683 0.545024 0.19393	0.768914 0.513774 0.380537	0.768914 0.513774 0.380537	0.815683 0.545024 0.19393	0.693681 0.693684 0.19393	0.768914 0.513774 0.380537	0.693681 0.693684 0.19393	0.653907 0.653909 0.380537	0.653907 0.653909 0.380537	0.693681 0.693684 0.19393	0.545021 0.815686 0.19393	0.653907 0.653909 0.380537	0.545021 0.815686 0.19393	0.513771 0.768915 0.380537	0.513771 0.768915 0.380537	0.545021 0.815686 0.19393	0.375417 0.906341 0.19393	0.513771 0.768915 0.380537	0.375417 0.906341 0.19393	0.353891 0.854373 0.380537	0.353891 0.854373 0.380537	0.375417 0.906341 0.19393	0.191385 0.962166 0.19393	0.353891 0.854373 0.380537	0.191385 0.962166 0.19393	0.180411 0.906997 0.380537	0.180411 0.906997 0.380537	0.191385 0.962166 0.19393	-7.6109e-007 0.981015 0.19393	0.180411 0.906997 0.380537	-7.6109e-007 0.981015 0.19393	-6.35304e-007 0.924766 0.380537	-7.6109e-007 0.981015 0.19393	-7.49269e-007 1.0 1.59009e-007	-0.19509 0.980785 1.46963e-007	-7.6109e-007 0.981015 0.19393	-0.19509 0.980785 1.46963e-007	-0.191387 0.962165 0.19393	-0.191387 0.962165 0.19393	-0.19509 0.980785 1.46963e-007	-0.382683 0.92388 1.80692e-007	-0.191387 0.962165 0.19393	-0.382683 0.92388 1.80692e-007	-0.375418 0.90634 0.19393	-0.375418 0.90634 0.19393	-0.382683 0.92388 1.80692e-007	-0.55557 0.83147 1.28893e-007	-0.375418 0.90634 0.19393	-0.55557 0.83147 1.28893e-007	-0.545023 0.815684 0.19393	-0.545023 0.815684 0.19393	-0.55557 0.83147 1.28893e-007	-0.707107 0.707107 1.96352e-007	-0.545023 0.815684 0.19393	-0.707107 0.707107 1.96352e-007	-0.693683 0.693683 0.19393	-0.693683 0.693683 0.19393	-0.707107 0.707107 1.96352e-007	-0.83147 0.55557 1.68646e-007	-0.693683 0.693683 0.19393	-0.83147 0.55557 1.68646e-007	-0.815684 0.545023 0.19393	-0.815684 0.545023 0.19393	-0.83147 0.55557 1.68646e-007	-0.92388 0.382683 1.66236e-007	-0.815684 0.545023 0.19393	-0.92388 0.382683 1.66236e-007	-0.90634 0.375418 0.19393	-0.90634 0.375418 0.19393	-0.92388 0.382683 1.66236e-007	-0.980785 0.19509 1.8551e-007	-0.90634 0.375418 0.19393	-0.980785 0.19509 1.8551e-007	-0.962165 0.191386 0.19393	-0.962165 0.191386 0.19393	-0.980785 0.19509 1.8551e-007	-1.0 -3.28859e-007 1.65032e-007	-0.962165 0.191386 0.19393	-1.0 -3.28859e-007 1.65032e-007	-0.981015 -3.15515e-007 0.19393	-0.981015 -3.15515e-007 0.19393	-1.0 -3.28859e-007 1.65032e-007	-0.980785 -0.195091 1.83101e-007	-0.981015 -3.15515e-007 0.19393	-0.980785 -0.195091 1.83101e-007	-0.962165 -0.191387 0.19393	-0.962165 -0.191387 0.19393	-0.980785 -0.195091 1.83101e-007	-0.923879 -0.382684 1.86715e-007	-0.962165 -0.191387 0.19393	-0.923879 -0.382684 1.86715e-007	-0.90634 -0.375419 0.19393	-0.90634 -0.375419 0.19393	-0.923879 -0.382684 1.86715e-007	-0.831469 -0.555571 2.27672e-007	-0.90634 -0.375419 0.19393	-0.831469 -0.555571 2.27672e-007	-0.815684 -0.545023 0.19393	-0.815684 -0.545023 0.19393	-0.831469 -0.555571 2.27672e-007	-0.707106 -0.707107 2.18035e-007	-0.815684 -0.545023 0.19393	-0.707106 -0.707107 2.18035e-007	-0.693682 -0.693683 0.19393	-0.693682 -0.693683 0.19393	-0.707106 -0.707107 2.18035e-007	-0.55557 -0.83147 1.78283e-007	-0.693682 -0.693683 0.19393	-0.55557 -0.83147 1.78283e-007	-0.545022 -0.815685 0.19393	-0.545022 -0.815685 0.19393	-0.55557 -0.83147 1.78283e-007	-0.382683 -0.92388 1.67441e-007	-0.545022 -0.815685 0.19393	-0.382683 -0.92388 1.67441e-007	-0.375418 -0.90634 0.19393	-0.375418 -0.90634 0.19393	-0.382683 -0.92388 1.67441e-007	-0.195089 -0.980785 1.67441e-007	-0.375418 -0.90634 0.19393	-0.195089 -0.980785 1.67441e-007	-0.191386 -0.962166 0.19393	-0.191386 -0.962166 0.19393	-0.195089 -0.980785 1.67441e-007	9.55258e-007 -1.0 1.59009e-007	-0.191386 -0.962166 0.19393	9.55258e-007 -1.0 1.59009e-007	9.51363e-007 -0.981015 0.19393	9.51363e-007 -0.981015 0.19393	9.55258e-007 -1.0 1.59009e-007	0.195091 -0.980785 1.79487e-007	9.51363e-007 -0.981015 0.19393	0.195091 -0.980785 1.79487e-007	0.191388 -0.962165 0.19393	0.191388 -0.962165 0.19393	0.195091 -0.980785 1.79487e-007	0.382684 -0.923879 1.79487e-007	0.191388 -0.962165 0.19393	0.382684 -0.923879 1.79487e-007	0.375419 -0.90634 0.19393	0.375419 -0.90634 0.19393	0.382684 -0.923879 1.79487e-007	0.555571 -0.831469 1.19257e-007	0.375419 -0.90634 0.19393	0.555571 -0.831469 1.19257e-007	0.545024 -0.815684 0.19393	0.545024 -0.815684 0.19393	0.555571 -0.831469 1.19257e-007	0.707108 -0.707106 0.0	0.545024 -0.815684 0.19393	0.707108 -0.707106 0.0	0.693684 -0.693682 0.19393	0.693684 -0.693682 0.19393	0.707108 -0.707106 0.0	0.83147 -0.555569 0.0	0.693684 -0.693682 0.19393	0.83147 -0.555569 0.0	0.815685 -0.545022 0.19393	0.815685 -0.545022 0.19393	0.83147 -0.555569 0.0	0.92388 -0.382682 1.44553e-007	0.815685 -0.545022 0.19393	0.92388 -0.382682 1.44553e-007	0.906341 -0.375417 0.19393	0.906341 -0.375417 0.19393	0.92388 -0.382682 1.44553e-007	0.980786 -0.195089 1.60213e-007	0.906341 -0.375417 0.19393	0.980786 -0.195089 1.60213e-007	0.962166 -0.191385 0.19393	0.962166 -0.191385 0.19393	0.980786 -0.195089 1.60213e-007	1.0 1.58768e-006 1.49372e-007	0.962166 -0.191385 0.19393	1.0 1.58768e-006 1.49372e-007	0.981015 1.57758e-006 0.19393	0.981015 1.57758e-006 0.19393	1.0 1.58768e-006 1.49372e-007	0.980785 0.195092 1.83101e-007	0.981015 1.57758e-006 0.19393	0.980785 0.195092 1.83101e-007	0.962165 0.191388 0.19393	0.962165 0.191388 0.19393	0.980785 0.195092 1.83101e-007	0.923879 0.382685 1.74669e-007	0.962165 0.191388 0.19393	0.923879 0.382685 1.74669e-007	0.906339 0.37542 0.19393	0.906339 0.37542 0.19393	0.923879 0.382685 1.74669e-007	0.831469 0.555572 1.60213e-007	0.906339 0.37542 0.19393	0.831469 0.555572 1.60213e-007	0.815683 0.545024 0.19393	0.815683 0.545024 0.19393	0.831469 0.555572 1.60213e-007	0.707105 0.707108 1.59009e-007	0.815683 0.545024 0.19393	0.707105 0.707108 1.59009e-007	0.693681 0.693684 0.19393	0.693681 0.693684 0.19393	0.707105 0.707108 1.59009e-007	0.555569 0.831471 2.57787e-007	0.693681 0.693684 0.19393	0.555569 0.831471 2.57787e-007	0.545021 0.815686 0.19393	0.545021 0.815686 0.19393	0.555569 0.831471 2.57787e-007	0.382682 0.92388 1.66236e-007	0.545021 0.815686 0.19393	0.382682 0.92388 1.66236e-007	0.375417 0.906341 0.19393	0.375417 0.906341 0.19393	0.382682 0.92388 1.66236e-007	0.195089 0.980786 1.31303e-007	0.375417 0.906341 0.19393	0.195089 0.980786 1.31303e-007	0.191385 0.962166 0.19393	0.191385 0.962166 0.19393	0.195089 0.980786 1.31303e-007	-7.49269e-007 1.0 1.59009e-007	0.191385 0.962166 0.19393	-7.49269e-007 1.0 1.59009e-007	-7.6109e-007 0.981015 0.19393	-7.49269e-007 1.0 1.59009e-007	-7.27371e-007 0.981015 -0.19393	-0.191387 0.962165 -0.19393	-7.49269e-007 1.0 1.59009e-007	-0.191387 0.962165 -0.19393	-0.19509 0.980785 1.46963e-007	-0.19509 0.980785 1.46963e-007	-0.191387 0.962165 -0.19393	-0.375418 0.90634 -0.19393	-0.19509 0.980785 1.46963e-007	-0.375418 0.90634 -0.19393	-0.382683 0.92388 1.80692e-007	-0.382683 0.92388 1.80692e-007	-0.375418 0.90634 -0.19393	-0.545023 0.815684 -0.19393	-0.382683 0.92388 1.80692e-007	-0.545023 0.815684 -0.19393	-0.55557 0.83147 1.28893e-007	-0.55557 0.83147 1.28893e-007	-0.545023 0.815684 -0.19393	-0.693683 0.693683 -0.19393	-0.55557 0.83147 1.28893e-007	-0.693683 0.693683 -0.19393	-0.707107 0.707107 1.96352e-007	-0.707107 0.707107 1.96352e-007	-0.693683 0.693683 -0.19393	-0.815684 0.545023 -0.19393	-0.707107 0.707107 1.96352e-007	-0.815684 0.545023 -0.19393	-0.83147 0.55557 1.68646e-007	-0.83147 0.55557 1.68646e-007	-0.815684 0.545023 -0.19393	-0.90634 0.375418 -0.19393	-0.83147 0.55557 1.68646e-007	-0.90634 0.375418 -0.19393	-0.92388 0.382683 1.66236e-007	-0.92388 0.382683 1.66236e-007	-0.90634 0.375418 -0.19393	-0.962165 0.191386 -0.19393	-0.92388 0.382683 1.66236e-007	-0.962165 0.191386 -0.19393	-0.980785 0.19509 1.8551e-007	-0.980785 0.19509 1.8551e-007	-0.962165 0.191386 -0.19393	-0.981015 -2.89022e-007 -0.19393	-0.980785 0.19509 1.8551e-007	-0.981015 -2.89022e-007 -0.19393	-1.0 -3.28859e-007 1.65032e-007	-1.0 -3.28859e-007 1.65032e-007	-0.981015 -2.89022e-007 -0.19393	-0.962165 -0.191387 -0.19393	-1.0 -3.28859e-007 1.65032e-007	-0.962165 -0.191387 -0.19393	-0.980785 -0.195091 1.83101e-007	-0.980785 -0.195091 1.83101e-007	-0.962165 -0.191387 -0.19393	-0.90634 -0.375419 -0.19393	-0.980785 -0.195091 1.83101e-007	-0.90634 -0.375419 -0.19393	-0.923879 -0.382684 1.86715e-007	-0.923879 -0.382684 1.86715e-007	-0.90634 -0.375419 -0.19393	-0.815684 -0.545023 -0.19393	-0.923879 -0.382684 1.86715e-007	-0.815684 -0.545023 -0.19393	-0.831469 -0.555571 2.27672e-007	-0.831469 -0.555571 2.27672e-007	-0.815684 -0.545023 -0.19393	-0.693682 -0.693683 -0.19393	-0.831469 -0.555571 2.27672e-007	-0.693682 -0.693683 -0.19393	-0.707106 -0.707107 2.18035e-007	-0.707106 -0.707107 2.18035e-007	-0.693682 -0.693683 -0.19393	-0.545022 -0.815685 -0.19393	-0.707106 -0.707107 2.18035e-007	-0.545022 -0.815685 -0.19393	-0.55557 -0.83147 1.78283e-007	-0.55557 -0.83147 1.78283e-007	-0.545022 -0.815685 -0.19393	-0.375418 -0.90634 -0.19393	-0.55557 -0.83147 1.78283e-007	-0.375418 -0.90634 -0.19393	-0.382683 -0.92388 1.67441e-007	-0.382683 -0.92388 1.67441e-007	-0.375418 -0.90634 -0.19393	-0.191386 -0.962166 -0.19393	-0.382683 -0.92388 1.67441e-007	-0.191386 -0.962166 -0.19393	-0.195089 -0.980785 1.67441e-007	-0.195089 -0.980785 1.67441e-007	-0.191386 -0.962166 -0.19393	9.32095e-007 -0.981015 -0.19393	-0.195089 -0.980785 1.67441e-007	9.32095e-007 -0.981015 -0.19393	9.55258e-007 -1.0 1.59009e-007	9.55258e-007 -1.0 1.59009e-007	9.32095e-007 -0.981015 -0.19393	0.191388 -0.962165 -0.19393	9.55258e-007 -1.0 1.59009e-007	0.191388 -0.962165 -0.19393	0.195091 -0.980785 1.79487e-007	0.195091 -0.980785 1.79487e-007	0.191388 -0.962165 -0.19393	0.375419 -0.90634 -0.19393	0.195091 -0.980785 1.79487e-007	0.375419 -0.90634 -0.19393	0.382684 -0.923879 1.79487e-007	0.382684 -0.923879 1.79487e-007	0.375419 -0.90634 -0.19393	0.545024 -0.815684 -0.19393	0.382684 -0.923879 1.79487e-007	0.545024 -0.815684 -0.19393	0.555571 -0.831469 1.19257e-007	0.555571 -0.831469 1.19257e-007	0.545024 -0.815684 -0.19393	0.693683 -0.693682 -0.19393	0.555571 -0.831469 1.19257e-007	0.693683 -0.693682 -0.19393	0.707108 -0.707106 0.0	0.707108 -0.707106 0.0	0.693683 -0.693682 -0.19393	0.815685 -0.545022 -0.19393	0.707108 -0.707106 0.0	0.815685 -0.545022 -0.19393	0.83147 -0.555569 0.0	0.83147 -0.555569 0.0	0.815685 -0.545022 -0.19393	0.906341 -0.375417 -0.19393	0.83147 -0.555569 0.0	0.906341 -0.375417 -0.19393	0.92388 -0.382682 1.44553e-007	0.92388 -0.382682 1.44553e-007	0.906341 -0.375417 -0.19393	0.962166 -0.191385 -0.19393	0.92388 -0.382682 1.44553e-007	0.962166 -0.191385 -0.19393	0.980786 -0.195089 1.60213e-007	0.980786 -0.195089 1.60213e-007	0.962166 -0.191385 -0.19393	0.981015 1.57035e-006 -0.19393	0.980786 -0.195089 1.60213e-007	0.981015 1.57035e-006 -0.19393	1.0 1.58768e-006 1.49372e-007	1.0 1.58768e-006 1.49372e-007	0.981015 1.57035e-006 -0.19393	0.962165 0.191388 -0.19393	1.0 1.58768e-006 1.49372e-007	0.962165 0.191388 -0.19393	0.980785 0.195092 1.83101e-007	0.980785 0.195092 1.83101e-007	0.962165 0.191388 -0.19393	0.906339 0.37542 -0.19393	0.980785 0.195092 1.83101e-007	0.906339 0.37542 -0.19393	0.923879 0.382685 1.74669e-007	0.923879 0.382685 1.74669e-007	0.906339 0.37542 -0.19393	0.815683 0.545024 -0.19393	0.923879 0.382685 1.74669e-007	0.815683 0.545024 -0.19393	0.831469 0.555572 1.60213e-007	0.831469 0.555572 1.60213e-007	0.815683 0.545024 -0.19393	0.693681 0.693684 -0.19393	0.831469 0.555572 1.60213e-007	0.693681 0.693684 -0.19393	0.707105 0.707108 1.59009e-007	0.707105 0.707108 1.59009e-007	0.693681 0.693684 -0.19393	0.545021 0.815686 -0.19393	0.707105 0.707108 1.59009e-007	0.545021 0.815686 -0.19393	0.555569 0.831471 2.57787e-007	0.555569 0.831471 2.57787e-007	0.545021 0.815686 -0.19393	0.375417 0.906341 -0.19393	0.555569 0.831471 2.57787e-007	0.375417 0.906341 -0.19393	0.382682 0.92388 1.66236e-007	0.382682 0.92388 1.66236e-007	0.375417 0.906341 -0.19393	0.191385 0.962166 -0.19393	0.382682 0.92388 1.66236e-007	0.191385 0.962166 -0.19393	0.195089 0.980786 1.31303e-007	0.195089 0.980786 1.31303e-007	0.191385 0.962166 -0.19393	-7.27371e-007 0.981015 -0.19393	0.195089 0.980786 1.31303e-007	-7.27371e-007 0.981015 -0.19393	-7.49269e-007 1.0 1.59009e-007	-7.27371e-007 0.981015 -0.19393	-6.42524e-007 0.924766 -0.380537	-0.180413 0.906997 -0.380537	-7.27371e-007 0.981015 -0.19393	-0.180413 0.906997 -0.380537	-0.191387 0.962165 -0.19393	-0.191387 0.962165 -0.19393	-0.180413 0.906997 -0.380537	-0.353893 0.854372 -0.380537	-0.191387 0.962165 -0.19393	-0.353893 0.854372 -0.380537	-0.375418 0.90634 -0.19393	-0.375418 0.90634 -0.19393	-0.353893 0.854372 -0.380537	-0.513772 0.768915 -0.380537	-0.375418 0.90634 -0.19393	-0.513772 0.768915 -0.380537	-0.545023 0.815684 -0.19393	-0.545023 0.815684 -0.19393	-0.513772 0.768915 -0.380537	-0.653908 0.653908 -0.380537	-0.545023 0.815684 -0.19393	-0.653908 0.653908 -0.380537	-0.693683 0.693683 -0.19393	-0.693683 0.693683 -0.19393	-0.653908 0.653908 -0.380537	-0.768915 0.513772 -0.380537	-0.693683 0.693683 -0.19393	-0.768915 0.513772 -0.380537	-0.815684 0.545023 -0.19393	-0.815684 0.545023 -0.19393	-0.768915 0.513772 -0.380537	-0.854372 0.353892 -0.380537	-0.815684 0.545023 -0.19393	-0.854372 0.353892 -0.380537	-0.90634 0.375418 -0.19393	-0.90634 0.375418 -0.19393	-0.854372 0.353892 -0.380537	-0.906997 0.180413 -0.380537	-0.90634 0.375418 -0.19393	-0.906997 0.180413 -0.380537	-0.962165 0.191386 -0.19393	-0.962165 0.191386 -0.19393	-0.906997 0.180413 -0.380537	-0.924766 -1.91313e-007 -0.380537	-0.962165 0.191386 -0.19393	-0.924766 -1.91313e-007 -0.380537	-0.981015 -2.89022e-007 -0.19393	-0.981015 -2.89022e-007 -0.19393	-0.924766 -1.91313e-007 -0.380537	-0.906996 -0.180413 -0.380537	-0.981015 -2.89022e-007 -0.19393	-0.906996 -0.180413 -0.380537	-0.962165 -0.191387 -0.19393	-0.962165 -0.191387 -0.19393	-0.906996 -0.180413 -0.380537	-0.854372 -0.353893 -0.380537	-0.962165 -0.191387 -0.19393	-0.854372 -0.353893 -0.380537	-0.90634 -0.375419 -0.19393	-0.90634 -0.375419 -0.19393	-0.854372 -0.353893 -0.380537	-0.768914 -0.513773 -0.380537	-0.90634 -0.375419 -0.19393	-0.768914 -0.513773 -0.380537	-0.815684 -0.545023 -0.19393	-0.815684 -0.545023 -0.19393	-0.768914 -0.513773 -0.380537	-0.653908 -0.653908 -0.380537	-0.815684 -0.545023 -0.19393	-0.653908 -0.653908 -0.380537	-0.693682 -0.693683 -0.19393	-0.693682 -0.693683 -0.19393	-0.653908 -0.653908 -0.380537	-0.513772 -0.768915 -0.380537	-0.693682 -0.693683 -0.19393	-0.513772 -0.768915 -0.380537	-0.545022 -0.815685 -0.19393	-0.545022 -0.815685 -0.19393	-0.513772 -0.768915 -0.380537	-0.353892 -0.854372 -0.380537	-0.545022 -0.815685 -0.19393	-0.353892 -0.854372 -0.380537	-0.375418 -0.90634 -0.19393	-0.375418 -0.90634 -0.19393	-0.353892 -0.854372 -0.380537	-0.180412 -0.906997 -0.380537	-0.375418 -0.90634 -0.19393	-0.180412 -0.906997 -0.380537	-0.191386 -0.962166 -0.19393	-0.191386 -0.962166 -0.19393	-0.180412 -0.906997 -0.380537	8.02553e-007 -0.924766 -0.380537	-0.191386 -0.962166 -0.19393	8.02553e-007 -0.924766 -0.380537	9.32095e-007 -0.981015 -0.19393	9.32095e-007 -0.981015 -0.19393	8.02553e-007 -0.924766 -0.380537	0.180414 -0.906996 -0.380537	9.32095e-007 -0.981015 -0.19393	0.180414 -0.906996 -0.380537	0.191388 -0.962165 -0.19393	0.191388 -0.962165 -0.19393	0.180414 -0.906996 -0.380537	0.353893 -0.854372 -0.380537	0.191388 -0.962165 -0.19393	0.353893 -0.854372 -0.380537	0.375419 -0.90634 -0.19393	0.375419 -0.90634 -0.19393	0.353893 -0.854372 -0.380537	0.513773 -0.768914 -0.380537	0.375419 -0.90634 -0.19393	0.513773 -0.768914 -0.380537	0.545024 -0.815684 -0.19393	0.545024 -0.815684 -0.19393	0.513773 -0.768914 -0.380537	0.653909 -0.653907 -0.380537	0.545024 -0.815684 -0.19393	0.653909 -0.653907 -0.380537	0.693683 -0.693682 -0.19393	0.693683 -0.693682 -0.19393	0.653909 -0.653907 -0.380537	0.768915 -0.513771 -0.380537	0.693683 -0.693682 -0.19393	0.768915 -0.513771 -0.380537	0.815685 -0.545022 -0.19393	0.815685 -0.545022 -0.19393	0.768915 -0.513771 -0.380537	0.854373 -0.353891 -0.380537	0.815685 -0.545022 -0.19393	0.854373 -0.353891 -0.380537	0.906341 -0.375417 -0.19393	0.906341 -0.375417 -0.19393	0.854373 -0.353891 -0.380537	0.906997 -0.180411 -0.380537	0.906341 -0.375417 -0.19393	0.906997 -0.180411 -0.380537	0.962166 -0.191385 -0.19393	0.962166 -0.191385 -0.19393	0.906997 -0.180411 -0.380537	0.924766 1.43786e-006 -0.380537	0.962166 -0.191385 -0.19393	0.924766 1.43786e-006 -0.380537	0.981015 1.57035e-006 -0.19393	0.981015 1.57035e-006 -0.19393	0.924766 1.43786e-006 -0.380537	0.906996 0.180414 -0.380537	0.981015 1.57035e-006 -0.19393	0.906996 0.180414 -0.380537	0.962165 0.191388 -0.19393	0.962165 0.191388 -0.19393	0.906996 0.180414 -0.380537	0.854371 0.353894 -0.380537	0.962165 0.191388 -0.19393	0.854371 0.353894 -0.380537	0.906339 0.37542 -0.19393	0.906339 0.37542 -0.19393	0.854371 0.353894 -0.380537	0.768914 0.513774 -0.380537	0.906339 0.37542 -0.19393	0.768914 0.513774 -0.380537	0.815683 0.545024 -0.19393	0.815683 0.545024 -0.19393	0.768914 0.513774 -0.380537	0.653907 0.653909 -0.380537	0.815683 0.545024 -0.19393	0.653907 0.653909 -0.380537	0.693681 0.693684 -0.19393	0.693681 0.693684 -0.19393	0.653907 0.653909 -0.380537	0.513771 0.768916 -0.380537	0.693681 0.693684 -0.19393	0.513771 0.768916 -0.380537	0.545021 0.815686 -0.19393	0.545021 0.815686 -0.19393	0.513771 0.768916 -0.380537	0.353891 0.854373 -0.380537	0.545021 0.815686 -0.19393	0.353891 0.854373 -0.380537	0.375417 0.906341 -0.19393	0.375417 0.906341 -0.19393	0.353891 0.854373 -0.380537	0.180411 0.906997 -0.380537	0.375417 0.906341 -0.19393	0.180411 0.906997 -0.380537	0.191385 0.962166 -0.19393	0.191385 0.962166 -0.19393	0.180411 0.906997 -0.380537	-6.42524e-007 0.924766 -0.380537	0.191385 0.962166 -0.19393	-6.42524e-007 0.924766 -0.380537	-7.27371e-007 0.981015 -0.19393	-6.42524e-007 0.924766 -0.380537	-5.41936e-007 0.83334 -0.552761	-0.162576 0.817327 -0.552761	-6.42524e-007 0.924766 -0.380537	-0.162576 0.817327 -0.552761	-0.180413 0.906997 -0.380537	-0.180413 0.906997 -0.380537	-0.162576 0.817327 -0.552761	-0.318905 0.769906 -0.552761	-0.180413 0.906997 -0.380537	-0.318905 0.769906 -0.552761	-0.353893 0.854372 -0.380537	-0.353893 0.854372 -0.380537	-0.318905 0.769906 -0.552761	-0.462979 0.692897 -0.552761	-0.353893 0.854372 -0.380537	-0.462979 0.692897 -0.552761	-0.513772 0.768915 -0.380537	-0.513772 0.768915 -0.380537	-0.462979 0.692897 -0.552761	-0.58926 0.58926 -0.552761	-0.513772 0.768915 -0.380537	-0.58926 0.58926 -0.552761	-0.653908 0.653908 -0.380537	-0.653908 0.653908 -0.380537	-0.58926 0.58926 -0.552761	-0.692897 0.462979 -0.552761	-0.653908 0.653908 -0.380537	-0.692897 0.462979 -0.552761	-0.768915 0.513772 -0.380537	-0.768915 0.513772 -0.380537	-0.692897 0.462979 -0.552761	-0.769906 0.318905 -0.552761	-0.768915 0.513772 -0.380537	-0.769906 0.318905 -0.552761	-0.854372 0.353892 -0.380537	-0.854372 0.353892 -0.380537	-0.769906 0.318905 -0.552761	-0.817327 0.162576 -0.552761	-0.854372 0.353892 -0.380537	-0.817327 0.162576 -0.552761	-0.906997 0.180413 -0.380537	-0.906997 0.180413 -0.380537	-0.817327 0.162576 -0.552761	-0.83334 -2.90795e-007 -0.552761	-0.906997 0.180413 -0.380537	-0.83334 -2.90795e-007 -0.552761	-0.924766 -1.91313e-007 -0.380537	-0.924766 -1.91313e-007 -0.380537	-0.83334 -2.90795e-007 -0.552761	-0.817327 -0.162577 -0.552761	-0.924766 -1.91313e-007 -0.380537	-0.817327 -0.162577 -0.552761	-0.906996 -0.180413 -0.380537	-0.906996 -0.180413 -0.380537	-0.817327 -0.162577 -0.552761	-0.769905 -0.318906 -0.552761	-0.906996 -0.180413 -0.380537	-0.769905 -0.318906 -0.552761	-0.854372 -0.353893 -0.380537	-0.854372 -0.353893 -0.380537	-0.769905 -0.318906 -0.552761	-0.692896 -0.462979 -0.552761	-0.854372 -0.353893 -0.380537	-0.692896 -0.462979 -0.552761	-0.768914 -0.513773 -0.380537	-0.768914 -0.513773 -0.380537	-0.692896 -0.462979 -0.552761	-0.58926 -0.589261 -0.552761	-0.768914 -0.513773 -0.380537	-0.58926 -0.589261 -0.552761	-0.653908 -0.653908 -0.380537	-0.653908 -0.653908 -0.380537	-0.58926 -0.589261 -0.552761	-0.462978 -0.692897 -0.552761	-0.653908 -0.653908 -0.380537	-0.462978 -0.692897 -0.552761	-0.513772 -0.768915 -0.380537	-0.513772 -0.768915 -0.380537	-0.462978 -0.692897 -0.552761	-0.318905 -0.769906 -0.552761	-0.513772 -0.768915 -0.380537	-0.318905 -0.769906 -0.552761	-0.353892 -0.854372 -0.380537	-0.353892 -0.854372 -0.380537	-0.318905 -0.769906 -0.552761	-0.162576 -0.817327 -0.552761	-0.353892 -0.854372 -0.380537	-0.162576 -0.817327 -0.552761	-0.180412 -0.906997 -0.380537	-0.180412 -0.906997 -0.380537	-0.162576 -0.817327 -0.552761	6.46478e-007 -0.83334 -0.552761	-0.180412 -0.906997 -0.380537	6.46478e-007 -0.83334 -0.552761	8.02553e-007 -0.924766 -0.380537	8.02553e-007 -0.924766 -0.380537	6.46478e-007 -0.83334 -0.552761	0.162577 -0.817327 -0.552761	8.02553e-007 -0.924766 -0.380537	0.162577 -0.817327 -0.552761	0.180414 -0.906996 -0.380537	0.180414 -0.906996 -0.380537	0.162577 -0.817327 -0.552761	0.318906 -0.769905 -0.552761	0.180414 -0.906996 -0.380537	0.318906 -0.769905 -0.552761	0.353893 -0.854372 -0.380537	0.353893 -0.854372 -0.380537	0.318906 -0.769905 -0.552761	0.46298 -0.692896 -0.552761	0.353893 -0.854372 -0.380537	0.46298 -0.692896 -0.552761	0.513773 -0.768914 -0.380537	0.513773 -0.768914 -0.380537	0.46298 -0.692896 -0.552761	0.589261 -0.589259 -0.552761	0.513773 -0.768914 -0.380537	0.589261 -0.589259 -0.552761	0.653909 -0.653907 -0.380537	0.653909 -0.653907 -0.380537	0.589261 -0.589259 -0.552761	0.692897 -0.462978 -0.552761	0.653909 -0.653907 -0.380537	0.692897 -0.462978 -0.552761	0.768915 -0.513771 -0.380537	0.768915 -0.513771 -0.380537	0.692897 -0.462978 -0.552761	0.769906 -0.318904 -0.552761	0.768915 -0.513771 -0.380537	0.769906 -0.318904 -0.552761	0.854373 -0.353891 -0.380537	0.854373 -0.353891 -0.380537	0.769906 -0.318904 -0.552761	0.817328 -0.162575 -0.552761	0.854373 -0.353891 -0.380537	0.817328 -0.162575 -0.552761	0.906997 -0.180411 -0.380537	0.906997 -0.180411 -0.380537	0.817328 -0.162575 -0.552761	0.83334 1.27974e-006 -0.552761	0.906997 -0.180411 -0.380537	0.83334 1.27974e-006 -0.552761	0.924766 1.43786e-006 -0.380537	0.924766 1.43786e-006 -0.380537	0.83334 1.27974e-006 -0.552761	0.817327 0.162578 -0.552761	0.924766 1.43786e-006 -0.380537	0.817327 0.162578 -0.552761	0.906996 0.180414 -0.380537	0.906996 0.180414 -0.380537	0.817327 0.162578 -0.552761	0.769905 0.318907 -0.552761	0.906996 0.180414 -0.380537	0.769905 0.318907 -0.552761	0.854371 0.353894 -0.380537	0.854371 0.353894 -0.380537	0.769905 0.318907 -0.552761	0.692896 0.46298 -0.552761	0.854371 0.353894 -0.380537	0.692896 0.46298 -0.552761	0.768914 0.513774 -0.380537	0.768914 0.513774 -0.380537	0.692896 0.46298 -0.552761	0.589259 0.589261 -0.552761	0.768914 0.513774 -0.380537	0.589259 0.589261 -0.552761	0.653907 0.653909 -0.380537	0.653907 0.653909 -0.380537	0.589259 0.589261 -0.552761	0.462977 0.692897 -0.552761	0.653907 0.653909 -0.380537	0.462977 0.692897 -0.552761	0.513771 0.768916 -0.380537	0.513771 0.768916 -0.380537	0.462977 0.692897 -0.552761	0.318904 0.769906 -0.552761	0.513771 0.768916 -0.380537	0.318904 0.769906 -0.552761	0.353891 0.854373 -0.380537	0.353891 0.854373 -0.380537	0.318904 0.769906 -0.552761	0.162575 0.817328 -0.552761	0.353891 0.854373 -0.380537	0.162575 0.817328 -0.552761	0.180411 0.906997 -0.380537	0.180411 0.906997 -0.380537	0.162575 0.817328 -0.552761	-5.41936e-007 0.83334 -0.552761	0.180411 0.906997 -0.380537	-5.41936e-007 0.83334 -0.552761	-6.42524e-007 0.924766 -0.380537	-5.41936e-007 0.83334 -0.552761	-4.18672e-007 0.710141 -0.704059	-0.138542 0.696496 -0.704059	-5.41936e-007 0.83334 -0.552761	-0.138542 0.696496 -0.704059	-0.162576 0.817327 -0.552761	-0.162576 0.817327 -0.552761	-0.138542 0.696496 -0.704059	-0.271759 0.656085 -0.704059	-0.162576 0.817327 -0.552761	-0.271759 0.656085 -0.704059	-0.318905 0.769906 -0.552761	-0.318905 0.769906 -0.552761	-0.271759 0.656085 -0.704059	-0.394533 0.590461 -0.704059	-0.318905 0.769906 -0.552761	-0.394533 0.590461 -0.704059	-0.462979 0.692897 -0.552761	-0.462979 0.692897 -0.552761	-0.394533 0.590461 -0.704059	-0.502146 0.502146 -0.704059	-0.462979 0.692897 -0.552761	-0.502146 0.502146 -0.704059	-0.58926 0.58926 -0.552761	-0.58926 0.58926 -0.552761	-0.502146 0.502146 -0.704059	-0.590461 0.394533 -0.704059	-0.58926 0.58926 -0.552761	-0.590461 0.394533 -0.704059	-0.692897 0.462979 -0.552761	-0.692897 0.462979 -0.552761	-0.590461 0.394533 -0.704059	-0.656085 0.271759 -0.704059	-0.692897 0.462979 -0.552761	-0.656085 0.271759 -0.704059	-0.769906 0.318905 -0.552761	-0.769906 0.318905 -0.552761	-0.656085 0.271759 -0.704059	-0.696496 0.138542 -0.704059	-0.769906 0.318905 -0.552761	-0.696496 0.138542 -0.704059	-0.817327 0.162576 -0.552761	-0.817327 0.162576 -0.552761	-0.696496 0.138542 -0.704059	-0.710141 -3.13105e-007 -0.704059	-0.817327 0.162576 -0.552761	-0.710141 -3.13105e-007 -0.704059	-0.83334 -2.90795e-007 -0.552761	-0.83334 -2.90795e-007 -0.552761	-0.710141 -3.13105e-007 -0.704059	-0.696496 -0.138542 -0.704059	-0.83334 -2.90795e-007 -0.552761	-0.696496 -0.138542 -0.704059	-0.817327 -0.162577 -0.552761	-0.817327 -0.162577 -0.552761	-0.696496 -0.138542 -0.704059	-0.656085 -0.27176 -0.704059	-0.817327 -0.162577 -0.552761	-0.656085 -0.27176 -0.704059	-0.769905 -0.318906 -0.552761	-0.769905 -0.318906 -0.552761	-0.656085 -0.27176 -0.704059	-0.590461 -0.394534 -0.704059	-0.769905 -0.318906 -0.552761	-0.590461 -0.394534 -0.704059	-0.692896 -0.462979 -0.552761	-0.692896 -0.462979 -0.552761	-0.590461 -0.394534 -0.704059	-0.502145 -0.502146 -0.704059	-0.692896 -0.462979 -0.552761	-0.502145 -0.502146 -0.704059	-0.58926 -0.589261 -0.552761	-0.58926 -0.589261 -0.552761	-0.502145 -0.502146 -0.704059	-0.394533 -0.590461 -0.704059	-0.58926 -0.589261 -0.552761	-0.394533 -0.590461 -0.704059	-0.462978 -0.692897 -0.552761	-0.462978 -0.692897 -0.552761	-0.394533 -0.590461 -0.704059	-0.271759 -0.656085 -0.704059	-0.462978 -0.692897 -0.552761	-0.271759 -0.656085 -0.704059	-0.318905 -0.769906 -0.552761	-0.318905 -0.769906 -0.552761	-0.271759 -0.656085 -0.704059	-0.138541 -0.696496 -0.704059	-0.318905 -0.769906 -0.552761	-0.138541 -0.696496 -0.704059	-0.162576 -0.817327 -0.552761	-0.162576 -0.817327 -0.552761	-0.138541 -0.696496 -0.704059	6.16612e-007 -0.710141 -0.704059	-0.162576 -0.817327 -0.552761	6.16612e-007 -0.710141 -0.704059	6.46478e-007 -0.83334 -0.552761	6.46478e-007 -0.83334 -0.552761	6.16612e-007 -0.710141 -0.704059	0.138542 -0.696496 -0.704059	6.46478e-007 -0.83334 -0.552761	0.138542 -0.696496 -0.704059	0.162577 -0.817327 -0.552761	0.162577 -0.817327 -0.552761	0.138542 -0.696496 -0.704059	0.27176 -0.656085 -0.704059	0.162577 -0.817327 -0.552761	0.27176 -0.656085 -0.704059	0.318906 -0.769905 -0.552761	0.318906 -0.769905 -0.552761	0.27176 -0.656085 -0.704059	0.394534 -0.59046 -0.704059	0.318906 -0.769905 -0.552761	0.394534 -0.59046 -0.704059	0.46298 -0.692896 -0.552761	0.46298 -0.692896 -0.552761	0.394534 -0.59046 -0.704059	0.502146 -0.502145 -0.704059	0.46298 -0.692896 -0.552761	0.502146 -0.502145 -0.704059	0.589261 -0.589259 -0.552761	0.589261 -0.589259 -0.552761	0.502146 -0.502145 -0.704059	0.590462 -0.394533 -0.704059	0.589261 -0.589259 -0.552761	0.590462 -0.394533 -0.704059	0.692897 -0.462978 -0.552761	0.692897 -0.462978 -0.552761	0.590462 -0.394533 -0.704059	0.656086 -0.271758 -0.704059	0.692897 -0.462978 -0.552761	0.656086 -0.271758 -0.704059	0.769906 -0.318904 -0.552761	0.769906 -0.318904 -0.552761	0.656086 -0.271758 -0.704059	0.696496 -0.138541 -0.704059	0.769906 -0.318904 -0.552761	0.696496 -0.138541 -0.704059	0.817328 -0.162575 -0.552761	0.817328 -0.162575 -0.552761	0.696496 -0.138541 -0.704059	0.710141 1.10126e-006 -0.704059	0.817328 -0.162575 -0.552761	0.710141 1.10126e-006 -0.704059	0.83334 1.27974e-006 -0.552761	0.83334 1.27974e-006 -0.552761	0.710141 1.10126e-006 -0.704059	0.696496 0.138543 -0.704059	0.83334 1.27974e-006 -0.552761	0.696496 0.138543 -0.704059	0.817327 0.162578 -0.552761	0.817327 0.162578 -0.552761	0.696496 0.138543 -0.704059	0.656085 0.27176 -0.704059	0.817327 0.162578 -0.552761	0.656085 0.27176 -0.704059	0.769905 0.318907 -0.552761	0.769905 0.318907 -0.552761	0.656085 0.27176 -0.704059	0.59046 0.394535 -0.704059	0.769905 0.318907 -0.552761	0.59046 0.394535 -0.704059	0.692896 0.46298 -0.552761	0.692896 0.46298 -0.552761	0.59046 0.394535 -0.704059	0.502145 0.502147 -0.704059	0.692896 0.46298 -0.552761	0.502145 0.502147 -0.704059	0.589259 0.589261 -0.552761	0.589259 0.589261 -0.552761	0.502145 0.502147 -0.704059	0.394532 0.590462 -0.704059	0.589259 0.589261 -0.552761	0.394532 0.590462 -0.704059	0.462977 0.692897 -0.552761	0.462977 0.692897 -0.552761	0.394532 0.590462 -0.704059	0.271758 0.656086 -0.704059	0.462977 0.692897 -0.552761	0.271758 0.656086 -0.704059	0.318904 0.769906 -0.552761	0.318904 0.769906 -0.552761	0.271758 0.656086 -0.704059	0.138541 0.696497 -0.704059	0.318904 0.769906 -0.552761	0.138541 0.696497 -0.704059	0.162575 0.817328 -0.552761	0.162575 0.817328 -0.552761	0.138541 0.696497 -0.704059	-4.18672e-007 0.710141 -0.704059	0.162575 0.817328 -0.552761	-4.18672e-007 0.710141 -0.704059	-5.41936e-007 0.83334 -0.552761	-4.18672e-007 0.710141 -0.704059	-3.32886e-007 0.559775 -0.828645	-0.109207 0.549019 -0.828645	-4.18672e-007 0.710141 -0.704059	-0.109207 0.549019 -0.828645	-0.138542 0.696496 -0.704059	-0.138542 0.696496 -0.704059	-0.109207 0.549019 -0.828645	-0.214217 0.517165 -0.828645	-0.138542 0.696496 -0.704059	-0.214217 0.517165 -0.828645	-0.271759 0.656085 -0.704059	-0.271759 0.656085 -0.704059	-0.214217 0.517165 -0.828645	-0.310994 0.465436 -0.828645	-0.271759 0.656085 -0.704059	-0.310994 0.465436 -0.828645	-0.394533 0.590461 -0.704059	-0.394533 0.590461 -0.704059	-0.310994 0.465436 -0.828645	-0.395821 0.395821 -0.828645	-0.394533 0.590461 -0.704059	-0.395821 0.395821 -0.828645	-0.502146 0.502146 -0.704059	-0.502146 0.502146 -0.704059	-0.395821 0.395821 -0.828645	-0.465436 0.310994 -0.828645	-0.502146 0.502146 -0.704059	-0.465436 0.310994 -0.828645	-0.590461 0.394533 -0.704059	-0.590461 0.394533 -0.704059	-0.465436 0.310994 -0.828645	-0.517165 0.214216 -0.828645	-0.590461 0.394533 -0.704059	-0.517165 0.214216 -0.828645	-0.656085 0.271759 -0.704059	-0.656085 0.271759 -0.704059	-0.517165 0.214216 -0.828645	-0.549019 0.109207 -0.828645	-0.656085 0.271759 -0.704059	-0.549019 0.109207 -0.828645	-0.696496 0.138542 -0.704059	-0.696496 0.138542 -0.704059	-0.549019 0.109207 -0.828645	-0.559775 -1.83806e-007 -0.828645	-0.696496 0.138542 -0.704059	-0.559775 -1.83806e-007 -0.828645	-0.710141 -3.13105e-007 -0.704059	-0.710141 -3.13105e-007 -0.704059	-0.559775 -1.83806e-007 -0.828645	-0.549019 -0.109207 -0.828645	-0.710141 -3.13105e-007 -0.704059	-0.549019 -0.109207 -0.828645	-0.696496 -0.138542 -0.704059	-0.696496 -0.138542 -0.704059	-0.549019 -0.109207 -0.828645	-0.517165 -0.214217 -0.828645	-0.696496 -0.138542 -0.704059	-0.517165 -0.214217 -0.828645	-0.656085 -0.27176 -0.704059	-0.656085 -0.27176 -0.704059	-0.517165 -0.214217 -0.828645	-0.465436 -0.310995 -0.828645	-0.656085 -0.27176 -0.704059	-0.465436 -0.310995 -0.828645	-0.590461 -0.394534 -0.704059	-0.590461 -0.394534 -0.704059	-0.465436 -0.310995 -0.828645	-0.39582 -0.395821 -0.828645	-0.590461 -0.394534 -0.704059	-0.39582 -0.395821 -0.828645	-0.502145 -0.502146 -0.704059	-0.502145 -0.502146 -0.704059	-0.39582 -0.395821 -0.828645	-0.310994 -0.465436 -0.828645	-0.502145 -0.502146 -0.704059	-0.310994 -0.465436 -0.828645	-0.394533 -0.590461 -0.704059	-0.394533 -0.590461 -0.704059	-0.310994 -0.465436 -0.828645	-0.214216 -0.517165 -0.828645	-0.394533 -0.590461 -0.704059	-0.214216 -0.517165 -0.828645	-0.271759 -0.656085 -0.704059	-0.271759 -0.656085 -0.704059	-0.214216 -0.517165 -0.828645	-0.109206 -0.549019 -0.828645	-0.271759 -0.656085 -0.704059	-0.109206 -0.549019 -0.828645	-0.138541 -0.696496 -0.704059	-0.138541 -0.696496 -0.704059	-0.109206 -0.549019 -0.828645	4.83762e-007 -0.559775 -0.828645	-0.138541 -0.696496 -0.704059	4.83762e-007 -0.559775 -0.828645	6.16612e-007 -0.710141 -0.704059	6.16612e-007 -0.710141 -0.704059	4.83762e-007 -0.559775 -0.828645	0.109207 -0.549019 -0.828645	6.16612e-007 -0.710141 -0.704059	0.109207 -0.549019 -0.828645	0.138542 -0.696496 -0.704059	0.138542 -0.696496 -0.704059	0.109207 -0.549019 -0.828645	0.214217 -0.517164 -0.828645	0.138542 -0.696496 -0.704059	0.214217 -0.517164 -0.828645	0.27176 -0.656085 -0.704059	0.27176 -0.656085 -0.704059	0.214217 -0.517164 -0.828645	0.310995 -0.465436 -0.828645	0.27176 -0.656085 -0.704059	0.310995 -0.465436 -0.828645	0.394534 -0.59046 -0.704059	0.394534 -0.59046 -0.704059	0.310995 -0.465436 -0.828645	0.395821 -0.39582 -0.828645	0.394534 -0.59046 -0.704059	0.395821 -0.39582 -0.828645	0.502146 -0.502145 -0.704059	0.502146 -0.502145 -0.704059	0.395821 -0.39582 -0.828645	0.465436 -0.310994 -0.828645	0.502146 -0.502145 -0.704059	0.465436 -0.310994 -0.828645	0.590462 -0.394533 -0.704059	0.590462 -0.394533 -0.704059	0.465436 -0.310994 -0.828645	0.517165 -0.214216 -0.828645	0.590462 -0.394533 -0.704059	0.517165 -0.214216 -0.828645	0.656086 -0.271758 -0.704059	0.656086 -0.271758 -0.704059	0.517165 -0.214216 -0.828645	0.549019 -0.109206 -0.828645	0.656086 -0.271758 -0.704059	0.549019 -0.109206 -0.828645	0.696496 -0.138541 -0.704059	0.696496 -0.138541 -0.704059	0.549019 -0.109206 -0.828645	0.559775 7.98088e-007 -0.828645	0.696496 -0.138541 -0.704059	0.559775 7.98088e-007 -0.828645	0.710141 1.10126e-006 -0.704059	0.710141 1.10126e-006 -0.704059	0.559775 7.98088e-007 -0.828645	0.549019 0.109208 -0.828645	0.710141 1.10126e-006 -0.704059	0.549019 0.109208 -0.828645	0.696496 0.138543 -0.704059	0.696496 0.138543 -0.704059	0.549019 0.109208 -0.828645	0.517164 0.214217 -0.828645	0.696496 0.138543 -0.704059	0.517164 0.214217 -0.828645	0.656085 0.27176 -0.704059	0.656085 0.27176 -0.704059	0.517164 0.214217 -0.828645	0.465435 0.310995 -0.828645	0.656085 0.27176 -0.704059	0.465435 0.310995 -0.828645	0.59046 0.394535 -0.704059	0.59046 0.394535 -0.704059	0.465435 0.310995 -0.828645	0.39582 0.395822 -0.828645	0.59046 0.394535 -0.704059	0.39582 0.395822 -0.828645	0.502145 0.502147 -0.704059	0.502145 0.502147 -0.704059	0.39582 0.395822 -0.828645	0.310993 0.465436 -0.828645	0.502145 0.502147 -0.704059	0.310993 0.465436 -0.828645	0.394532 0.590462 -0.704059	0.394532 0.590462 -0.704059	0.310993 0.465436 -0.828645	0.214216 0.517165 -0.828645	0.394532 0.590462 -0.704059	0.214216 0.517165 -0.828645	0.271758 0.656086 -0.704059	0.271758 0.656086 -0.704059	0.214216 0.517165 -0.828645	0.109206 0.549019 -0.828645	0.271758 0.656086 -0.704059	0.109206 0.549019 -0.828645	0.138541 0.696497 -0.704059	0.138541 0.696497 -0.704059	0.109206 0.549019 -0.828645	-3.32886e-007 0.559775 -0.828645	0.138541 0.696497 -0.704059	-3.32886e-007 0.559775 -0.828645	-4.18672e-007 0.710141 -0.704059	-3.32886e-007 0.559775 -0.828645	-2.25899e-007 0.387887 -0.921707	-0.0756729 0.380434 -0.921707	-3.32886e-007 0.559775 -0.828645	-0.0756729 0.380434 -0.921707	-0.109207 0.549019 -0.828645	-0.109207 0.549019 -0.828645	-0.0756729 0.380434 -0.921707	-0.148438 0.358361 -0.921707	-0.109207 0.549019 -0.828645	-0.148438 0.358361 -0.921707	-0.214217 0.517165 -0.828645	-0.214217 0.517165 -0.828645	-0.148438 0.358361 -0.921707	-0.215498 0.322516 -0.921707	-0.214217 0.517165 -0.828645	-0.215498 0.322516 -0.921707	-0.310994 0.465436 -0.828645	-0.310994 0.465436 -0.828645	-0.215498 0.322516 -0.921707	-0.274277 0.274277 -0.921707	-0.310994 0.465436 -0.828645	-0.274277 0.274277 -0.921707	-0.395821 0.395821 -0.828645	-0.395821 0.395821 -0.828645	-0.274277 0.274277 -0.921707	-0.322516 0.215498 -0.921707	-0.395821 0.395821 -0.828645	-0.322516 0.215498 -0.921707	-0.465436 0.310994 -0.828645	-0.465436 0.310994 -0.828645	-0.322516 0.215498 -0.921707	-0.358361 0.148438 -0.921707	-0.465436 0.310994 -0.828645	-0.358361 0.148438 -0.921707	-0.517165 0.214216 -0.828645	-0.517165 0.214216 -0.828645	-0.358361 0.148438 -0.921707	-0.380434 0.0756728 -0.921707	-0.517165 0.214216 -0.828645	-0.380434 0.0756728 -0.921707	-0.549019 0.109207 -0.828645	-0.549019 0.109207 -0.828645	-0.380434 0.0756728 -0.921707	-0.387887 0.0 -0.921707	-0.549019 0.109207 -0.828645	-0.387887 0.0 -0.921707	-0.559775 -1.83806e-007 -0.828645	-0.559775 -1.83806e-007 -0.828645	-0.387887 0.0 -0.921707	-0.380434 -0.0756731 -0.921707	-0.559775 -1.83806e-007 -0.828645	-0.380434 -0.0756731 -0.921707	-0.549019 -0.109207 -0.828645	-0.549019 -0.109207 -0.828645	-0.380434 -0.0756731 -0.921707	-0.358361 -0.148438 -0.921707	-0.549019 -0.109207 -0.828645	-0.358361 -0.148438 -0.921707	-0.517165 -0.214217 -0.828645	-0.517165 -0.214217 -0.828645	-0.358361 -0.148438 -0.921707	-0.322516 -0.215498 -0.921707	-0.517165 -0.214217 -0.828645	-0.322516 -0.215498 -0.921707	-0.465436 -0.310995 -0.828645	-0.465436 -0.310995 -0.828645	-0.322516 -0.215498 -0.921707	-0.274277 -0.274278 -0.921707	-0.465436 -0.310995 -0.828645	-0.274277 -0.274278 -0.921707	-0.39582 -0.395821 -0.828645	-0.39582 -0.395821 -0.828645	-0.274277 -0.274278 -0.921707	-0.215498 -0.322516 -0.921707	-0.39582 -0.395821 -0.828645	-0.215498 -0.322516 -0.921707	-0.310994 -0.465436 -0.828645	-0.310994 -0.465436 -0.828645	-0.215498 -0.322516 -0.921707	-0.148438 -0.358361 -0.921707	-0.310994 -0.465436 -0.828645	-0.148438 -0.358361 -0.921707	-0.214216 -0.517165 -0.828645	-0.214216 -0.517165 -0.828645	-0.148438 -0.358361 -0.921707	-0.0756726 -0.380434 -0.921707	-0.214216 -0.517165 -0.828645	-0.0756726 -0.380434 -0.921707	-0.109206 -0.549019 -0.828645	-0.109206 -0.549019 -0.828645	-0.0756726 -0.380434 -0.921707	3.90841e-007 -0.387887 -0.921707	-0.109206 -0.549019 -0.828645	3.90841e-007 -0.387887 -0.921707	4.83762e-007 -0.559775 -0.828645	4.83762e-007 -0.559775 -0.828645	3.90841e-007 -0.387887 -0.921707	0.0756733 -0.380434 -0.921707	4.83762e-007 -0.559775 -0.828645	0.0756733 -0.380434 -0.921707	0.109207 -0.549019 -0.828645	0.109207 -0.549019 -0.828645	0.0756733 -0.380434 -0.921707	0.148438 -0.35836 -0.921707	0.109207 -0.549019 -0.828645	0.148438 -0.35836 -0.921707	0.214217 -0.517164 -0.828645	0.214217 -0.517164 -0.828645	0.148438 -0.35836 -0.921707	0.215499 -0.322516 -0.921707	0.214217 -0.517164 -0.828645	0.215499 -0.322516 -0.921707	0.310995 -0.465436 -0.828645	0.310995 -0.465436 -0.828645	0.215499 -0.322516 -0.921707	0.274278 -0.274277 -0.921707	0.310995 -0.465436 -0.828645	0.274278 -0.274277 -0.921707	0.395821 -0.39582 -0.828645	0.395821 -0.39582 -0.828645	0.274278 -0.274277 -0.921707	0.322516 -0.215498 -0.921707	0.395821 -0.39582 -0.828645	0.322516 -0.215498 -0.921707	0.465436 -0.310994 -0.828645	0.465436 -0.310994 -0.828645	0.322516 -0.215498 -0.921707	0.358361 -0.148437 -0.921707	0.465436 -0.310994 -0.828645	0.358361 -0.148437 -0.921707	0.517165 -0.214216 -0.828645	0.517165 -0.214216 -0.828645	0.358361 -0.148437 -0.921707	0.380434 -0.0756724 -0.921707	0.517165 -0.214216 -0.828645	0.380434 -0.0756724 -0.921707	0.549019 -0.109206 -0.828645	0.549019 -0.109206 -0.828645	0.380434 -0.0756724 -0.921707	0.387887 5.52197e-007 -0.921707	0.549019 -0.109206 -0.828645	0.387887 5.52197e-007 -0.921707	0.559775 7.98088e-007 -0.828645	0.559775 7.98088e-007 -0.828645	0.387887 5.52197e-007 -0.921707	0.380433 0.0756736 -0.921707	0.559775 7.98088e-007 -0.828645	0.380433 0.0756736 -0.921707	0.549019 0.109208 -0.828645	0.549019 0.109208 -0.828645	0.380433 0.0756736 -0.921707	0.35836 0.148438 -0.921707	0.549019 0.109208 -0.828645	0.35836 0.148438 -0.921707	0.517164 0.214217 -0.828645	0.517164 0.214217 -0.828645	0.35836 0.148438 -0.921707	0.322516 0.215499 -0.921707	0.517164 0.214217 -0.828645	0.322516 0.215499 -0.921707	0.465435 0.310995 -0.828645	0.465435 0.310995 -0.828645	0.322516 0.215499 -0.921707	0.274277 0.274278 -0.921707	0.465435 0.310995 -0.828645	0.274277 0.274278 -0.921707	0.39582 0.395822 -0.828645	0.39582 0.395822 -0.828645	0.274277 0.274278 -0.921707	0.215498 0.322517 -0.921707	0.39582 0.395822 -0.828645	0.215498 0.322517 -0.921707	0.310993 0.465436 -0.828645	0.310993 0.465436 -0.828645	0.215498 0.322517 -0.921707	0.148437 0.358361 -0.921707	0.310993 0.465436 -0.828645	0.148437 0.358361 -0.921707	0.214216 0.517165 -0.828645	0.214216 0.517165 -0.828645	0.148437 0.358361 -0.921707	0.0756723 0.380434 -0.921707	0.214216 0.517165 -0.828645	0.0756723 0.380434 -0.921707	0.109206 0.549019 -0.828645	0.109206 0.549019 -0.828645	0.0756723 0.380434 -0.921707	-2.25899e-007 0.387887 -0.921707	0.109206 0.549019 -0.828645	-2.25899e-007 0.387887 -0.921707	-3.32886e-007 0.559775 -0.828645	-2.25899e-007 0.387887 -0.921707	-1.48109e-007 0.200969 -0.979598	-0.039207 0.197107 -0.979598	-2.25899e-007 0.387887 -0.921707	-0.039207 0.197107 -0.979598	-0.0756729 0.380434 -0.921707	-0.0756729 0.380434 -0.921707	-0.039207 0.197107 -0.979598	-0.0769073 0.185671 -0.979598	-0.0756729 0.380434 -0.921707	-0.0769073 0.185671 -0.979598	-0.148438 0.358361 -0.921707	-0.148438 0.358361 -0.921707	-0.0769073 0.185671 -0.979598	-0.111652 0.167099 -0.979598	-0.148438 0.358361 -0.921707	-0.111652 0.167099 -0.979598	-0.215498 0.322516 -0.921707	-0.215498 0.322516 -0.921707	-0.111652 0.167099 -0.979598	-0.142106 0.142106 -0.979598	-0.215498 0.322516 -0.921707	-0.142106 0.142106 -0.979598	-0.274277 0.274277 -0.921707	-0.274277 0.274277 -0.921707	-0.142106 0.142106 -0.979598	-0.167099 0.111652 -0.979598	-0.274277 0.274277 -0.921707	-0.167099 0.111652 -0.979598	-0.322516 0.215498 -0.921707	-0.322516 0.215498 -0.921707	-0.167099 0.111652 -0.979598	-0.185671 0.0769073 -0.979598	-0.322516 0.215498 -0.921707	-0.185671 0.0769073 -0.979598	-0.358361 0.148438 -0.921707	-0.358361 0.148438 -0.921707	-0.185671 0.0769073 -0.979598	-0.197107 0.039207 -0.979598	-0.358361 0.148438 -0.921707	-0.197107 0.039207 -0.979598	-0.380434 0.0756728 -0.921707	-0.380434 0.0756728 -0.921707	-0.197107 0.039207 -0.979598	-0.200969 0.0 -0.979598	-0.380434 0.0756728 -0.921707	-0.200969 0.0 -0.979598	-0.387887 0.0 -0.921707	-0.387887 0.0 -0.921707	-0.200969 0.0 -0.979598	-0.197107 -0.0392071 -0.979598	-0.387887 0.0 -0.921707	-0.197107 -0.0392071 -0.979598	-0.380434 -0.0756731 -0.921707	-0.380434 -0.0756731 -0.921707	-0.197107 -0.0392071 -0.979598	-0.185671 -0.0769074 -0.979598	-0.380434 -0.0756731 -0.921707	-0.185671 -0.0769074 -0.979598	-0.358361 -0.148438 -0.921707	-0.358361 -0.148438 -0.921707	-0.185671 -0.0769074 -0.979598	-0.167099 -0.111652 -0.979598	-0.358361 -0.148438 -0.921707	-0.167099 -0.111652 -0.979598	-0.322516 -0.215498 -0.921707	-0.322516 -0.215498 -0.921707	-0.167099 -0.111652 -0.979598	-0.142106 -0.142106 -0.979598	-0.322516 -0.215498 -0.921707	-0.142106 -0.142106 -0.979598	-0.274277 -0.274278 -0.921707	-0.274277 -0.274278 -0.921707	-0.142106 -0.142106 -0.979598	-0.111652 -0.167099 -0.979598	-0.274277 -0.274278 -0.921707	-0.111652 -0.167099 -0.979598	-0.215498 -0.322516 -0.921707	-0.215498 -0.322516 -0.921707	-0.111652 -0.167099 -0.979598	-0.0769072 -0.185671 -0.979598	-0.215498 -0.322516 -0.921707	-0.0769072 -0.185671 -0.979598	-0.148438 -0.358361 -0.921707	-0.148438 -0.358361 -0.921707	-0.0769072 -0.185671 -0.979598	-0.0392068 -0.197107 -0.979598	-0.148438 -0.358361 -0.921707	-0.0392068 -0.197107 -0.979598	-0.0756726 -0.380434 -0.921707	-0.0756726 -0.380434 -0.921707	-0.0392068 -0.197107 -0.979598	2.07622e-007 -0.200969 -0.979598	-0.0756726 -0.380434 -0.921707	2.07622e-007 -0.200969 -0.979598	3.90841e-007 -0.387887 -0.921707	3.90841e-007 -0.387887 -0.921707	2.07622e-007 -0.200969 -0.979598	0.0392072 -0.197107 -0.979598	3.90841e-007 -0.387887 -0.921707	0.0392072 -0.197107 -0.979598	0.0756733 -0.380434 -0.921707	0.0756733 -0.380434 -0.921707	0.0392072 -0.197107 -0.979598	0.0769076 -0.185671 -0.979598	0.0756733 -0.380434 -0.921707	0.0769076 -0.185671 -0.979598	0.148438 -0.35836 -0.921707	0.148438 -0.35836 -0.921707	0.0769076 -0.185671 -0.979598	0.111652 -0.167099 -0.979598	0.148438 -0.35836 -0.921707	0.111652 -0.167099 -0.979598	0.215499 -0.322516 -0.921707	0.215499 -0.322516 -0.921707	0.111652 -0.167099 -0.979598	0.142106 -0.142106 -0.979598	0.215499 -0.322516 -0.921707	0.142106 -0.142106 -0.979598	0.274278 -0.274277 -0.921707	0.274278 -0.274277 -0.921707	0.142106 -0.142106 -0.979598	0.167099 -0.111652 -0.979598	0.274278 -0.274277 -0.921707	0.167099 -0.111652 -0.979598	0.322516 -0.215498 -0.921707	0.322516 -0.215498 -0.921707	0.167099 -0.111652 -0.979598	0.185671 -0.0769071 -0.979598	0.322516 -0.215498 -0.921707	0.185671 -0.0769071 -0.979598	0.358361 -0.148437 -0.921707	0.358361 -0.148437 -0.921707	0.185671 -0.0769071 -0.979598	0.197107 -0.0392067 -0.979598	0.358361 -0.148437 -0.921707	0.197107 -0.0392067 -0.979598	0.380434 -0.0756724 -0.921707	0.380434 -0.0756724 -0.921707	0.197107 -0.0392067 -0.979598	0.200969 2.84734e-007 -0.979598	0.380434 -0.0756724 -0.921707	0.200969 2.84734e-007 -0.979598	0.387887 5.52197e-007 -0.921707	0.387887 5.52197e-007 -0.921707	0.200969 2.84734e-007 -0.979598	0.197107 0.0392074 -0.979598	0.387887 5.52197e-007 -0.921707	0.197107 0.0392074 -0.979598	0.380433 0.0756736 -0.921707	0.380433 0.0756736 -0.921707	0.197107 0.0392074 -0.979598	0.185671 0.0769077 -0.979598	0.380433 0.0756736 -0.921707	0.185671 0.0769077 -0.979598	0.35836 0.148438 -0.921707	0.35836 0.148438 -0.921707	0.185671 0.0769077 -0.979598	0.167099 0.111652 -0.979598	0.35836 0.148438 -0.921707	0.167099 0.111652 -0.979598	0.322516 0.215499 -0.921707	0.322516 0.215499 -0.921707	0.167099 0.111652 -0.979598	0.142106 0.142107 -0.979598	0.322516 0.215499 -0.921707	0.142106 0.142107 -0.979598	0.274277 0.274278 -0.921707	0.274277 0.274278 -0.921707	0.142106 0.142107 -0.979598	0.111652 0.167099 -0.979598	0.274277 0.274278 -0.921707	0.111652 0.167099 -0.979598	0.215498 0.322517 -0.921707	0.215498 0.322517 -0.921707	0.111652 0.167099 -0.979598	0.076907 0.185671 -0.979598	0.215498 0.322517 -0.921707	0.076907 0.185671 -0.979598	0.148437 0.358361 -0.921707	0.148437 0.358361 -0.921707	0.076907 0.185671 -0.979598	0.0392067 0.197107 -0.979598	0.148437 0.358361 -0.921707	0.0392067 0.197107 -0.979598	0.0756723 0.380434 -0.921707	0.0756723 0.380434 -0.921707	0.0392067 0.197107 -0.979598	-1.48109e-007 0.200969 -0.979598	0.0756723 0.380434 -0.921707	-1.48109e-007 0.200969 -0.979598	-2.25899e-007 0.387887 -0.921707	0.0 0.0 -1.0	-0.039207 0.197107 -0.979598	-1.48109e-007 0.200969 -0.979598	0.0 0.0 -1.0	-0.0769073 0.185671 -0.979598	-0.039207 0.197107 -0.979598	0.0 0.0 -1.0	-0.111652 0.167099 -0.979598	-0.0769073 0.185671 -0.979598	0.0 0.0 -1.0	-0.142106 0.142106 -0.979598	-0.111652 0.167099 -0.979598	0.0 0.0 -1.0	-0.167099 0.111652 -0.979598	-0.142106 0.142106 -0.979598	0.0 0.0 -1.0	-0.185671 0.0769073 -0.979598	-0.167099 0.111652 -0.979598	0.0 0.0 -1.0	-0.197107 0.039207 -0.979598	-0.185671 0.0769073 -0.979598	0.0 0.0 -1.0	-0.200969 0.0 -0.979598	-0.197107 0.039207 -0.979598	0.0 0.0 -1.0	-0.197107 -0.0392071 -0.979598	-0.200969 0.0 -0.979598	0.0 0.0 -1.0	-0.185671 -0.0769074 -0.979598	-0.197107 -0.0392071 -0.979598	0.0 0.0 -1.0	-0.167099 -0.111652 -0.979598	-0.185671 -0.0769074 -0.979598	0.0 0.0 -1.0	-0.142106 -0.142106 -0.979598	-0.167099 -0.111652 -0.979598	0.0 0.0 -1.0	-0.111652 -0.167099 -0.979598	-0.142106 -0.142106 -0.979598	0.0 0.0 -1.0	-0.0769072 -0.185671 -0.979598	-0.111652 -0.167099 -0.979598	0.0 0.0 -1.0	-0.0392068 -0.197107 -0.979598	-0.0769072 -0.185671 -0.979598	0.0 0.0 -1.0	2.07622e-007 -0.200969 -0.979598	-0.0392068 -0.197107 -0.979598	0.0 0.0 -1.0	0.0392072 -0.197107 -0.979598	2.07622e-007 -0.200969 -0.979598	0.0 0.0 -1.0	0.0769076 -0.185671 -0.979598	0.0392072 -0.197107 -0.979598	0.0 0.0 -1.0	0.111652 -0.167099 -0.979598	0.0769076 -0.185671 -0.979598	0.0 0.0 -1.0	0.142106 -0.142106 -0.979598	0.111652 -0.167099 -0.979598	0.0 0.0 -1.0	0.167099 -0.111652 -0.979598	0.142106 -0.142106 -0.979598	0.0 0.0 -1.0	0.185671 -0.0769071 -0.979598	0.167099 -0.111652 -0.979598	0.0 0.0 -1.0	0.197107 -0.0392067 -0.979598	0.185671 -0.0769071 -0.979598	0.0 0.0 -1.0	0.200969 2.84734e-007 -0.979598	0.197107 -0.0392067 -0.979598	0.0 0.0 -1.0	0.197107 0.0392074 -0.979598	0.200969 2.84734e-007 -0.979598	0.0 0.0 -1.0	0.185671 0.0769077 -0.979598	0.197107 0.0392074 -0.979598	0.0 0.0 -1.0	0.167099 0.111652 -0.979598	0.185671 0.0769077 -0.979598	0.0 0.0 -1.0	0.142106 0.142107 -0.979598	0.167099 0.111652 -0.979598	0.0 0.0 -1.0	0.111652 0.167099 -0.979598	0.142106 0.142107 -0.979598	0.0 0.0 -1.0	0.076907 0.185671 -0.979598	0.111652 0.167099 -0.979598	0.0 0.0 -1.0	0.0392067 0.197107 -0.979598	0.076907 0.185671 -0.979598	0.0 0.0 -1.0	-1.48109e-007 0.200969 -0.979598	0.0392067 0.197107 -0.979598]\n" to:out_file
		format "AttributeEnd \n" to:out_file
		--**** End Object: Sphere08 **** 

			
		format "WorldEnd \n" to:out_file

		
		
			
			setWaitCursor()
			theCommand = "\"\"" + PBRT_EXE_PATH + "pbrt.exe\" \"" + fname  + ".pbrt\"\""			

			HiddenDOSCommand theCommand startpath:SAVEPATH --donotwait:true

			setArrowCursor()	

			close out_file
			
			-- display the image
--	theBitmap = undefined
			theBitmap = (openBitMap (SAVEPATH + fname + ".exr"))
				
--			dispBmp = copy theBitmap
	--		close theBitmap

				
			display theBitmap

	--		close theBitmap
	
	)
		
--// Matte
plugin material PBRT_Matte name:"PBRT Matte Mtl" classID:#(0x1b0818ca, 0x1ad4c418) extends:Standard replaceUI:true version:1
( 
	Local SigmaTmp = 0.0
	Local SigmaOff	= 0.0

	parameters main rollout:params
	(
		Kd type:#color default:[140,140,140] ui:matte_Diffuse_col 
		KdTexture type:#textureMap ui:matte_Texture_btn
		KdTextureEnable type:#boolean default:false ui:matte_Tex_chk
		
		alphaAmount type:#float default:1.0 ui:alpha_spn
		alphaEnable type:#boolean default:false ui:alphaEnabled_chk
		alphaImage type:#textureMap ui:alphaMap_btn
		
		sigmaMap type:#textureMap ui:matte_Sigma_Tex	
		sigma type:#float default:0.0 ui:matte_Sigma_spn
		sigmaMapEnabled type:#boolean default:false ui:matte_Sigma_chk
		
		mode type:#radiobtnIndex default:1 ui:matte_Mode 
		
		bumpMap type:#textureMap ui:matte_Bumpmap
		bumpMapEnable type:#boolean default:false ui:matte_Bump_chk
		bumpMapAmount type:#float default:0.05 ui:matte_Bump_spn
		
		LightColorEnable type:#boolean default:false ui:LightSamplesEnabled_chk
		LightColor type:#color default:[255,255,255] ui:lightColor_col
--		LightCMap type:#textureMap ui:LightSamplesTexture_btn
		LightMult type:#float default:1.0 ui:LightMult_spn
		LightSamples type:#integer default:1 ui:LightSamples_spn
		
		on Kd set val do delegate.diffuse = val
		on KdTexture set val do delegate.diffuseMap = val	
		on KdTextureEnable set val do delegate.diffuseMapEnable = val
			
		on alphaAmount set val do delegate.opacityMapAmount = val	
		on alphaEnable set val do delegate.opacityMapEnable = val	
		on alphaImage set val do delegate.opacityMap = val		
			
		on bumpMap set val do delegate.bumpMap = val
		on bumpMapAmount set val do delegate.bumpMapAmount = val
		on bumpMapEnable set val do delegate.bumpMapEnable = val	
			
		on LightColorEnable set val do delegate.useSelfIllumColor = val		
		on LightColor set val do delegate.selfIllumColor = val
--		on LightCMap set val do delegate.selfillumMap = val		
		on LightMult set val do delegate.selfillumMapAmount = val
	)

	
	rollout params "Matte Parameters"
	(
		GroupBox matte_group "" pos:[5,0] width:315 height:98 
		colorPicker matte_Diffuse_col "Diffuse Color: " pos:[15,13] width:110 height:20
		button matte_GetMapColors "Update Color from Bitmap" pos:[150,15] width:130 height:17 enabled:false toolTip:"Get Diffuse Color from Bitmap"
		label matte_dif_tex_lbl "Diffuse Texture:" pos:[14,48] width:100 height:20
		checkbox matte_Tex_chk "" pos:[98,48] width:15 height:15 checked:false	
		mapButton matte_Texture_btn "None" pos:[117,46] width:195 height:20

		spinner 		alpha_spn "Alpha:" pos:[34,75] width:58 height:16 range:[0,1,1] type:#integer
		checkbox 		alphaEnabled_chk "" pos:[100,75]
		mapButton 	alphaMap_btn "Alpah Map" pos:[120,73] width:195 height:20
		
		GroupBox matte_group2 "Mode:" pos:[5,102] width:315 height:70	
		radiobuttons matte_Mode "" pos:[70,119] width:193 height:16 labels:#("Lambertian", "Oren-Nayar") columns:2
		spinner matte_Sigma_spn "Sigma:" pos:[35,146] width:70 height:16 enabled:false range:[0,100,0]
		checkbox matte_Sigma_chk "" pos:[116,146] width:15 height:15 checked:false
		mapButton matte_Sigma_Tex "None" pos:[140,144] width:175 height:20 enabled:false

		GroupBox matte_group3 "" pos:[5,174] width:315 height:40			
		spinner matte_Bump_spn "Bumpmap:" pos:[48,189] width:65 height:16 range:[-1,1,0.05]
		checkbox matte_Bump_chk "" pos:[116,189] width:15 height:15 checked:false	
		mapButton matte_Bumpmap "T" pos:[135,187] width:180 height:20
		
		GroupBox LightSamples_gup "Emitter" pos:[6,218] width:315 height:84
		checkbox LightSamplesEnabled_chk "On/Off:" pos:[18,234] width:57 height:13 	
		colorPicker lightColor_col "Light Color" pos:[19,250] width:92 height:20 enabled:false
		spinner LightMult_spn "Light Multiplier:" pos:[70,280] width:80 height:16 enabled:false range:[0,1000,1]
		spinner LightSamples_spn "Light Samples:" pos:[210,281] width:80 height:16 enabled:false range:[0,1000,4] type:#integer	

		
		GroupBox grp_Preview "" pos:[6,304] width:315 height:75				
		
		radiobuttons Preview_Obj "" pos:[15,317] width:209 height:16 labels:#("Box", "Sphere", "Torus", "Plane", "Surface") default:2 columns:5
		radiobuttons Preview_Type "" pos:[100,340] width:209 height:16 labels:#("Direct Lighting", "Metropolis") default:1 columns:2 
		radiobuttons Preview_Quality "" pos:[100,360] width:153 height:16 labels:#("Low", "Med", "High", "Best") default:2 columns:4
		button preview "Preview" pos:[15,340] width:62 height:35
	
		on params open do (
			if KdTexture != undefined do (
				matte_GetMapColors.enabled = true
				
					-- try to update the buttons text
					try(matte_Texture_btn.text = (filenameFromPath KdTexture.PBRT_BitMap.filename)
						matte_Texture_btn.tooltip = KdTexture.PBRT_BitMap.filename)catch()
			)
			
			if matte_Mode.state == 2 then (
				
					matte_Sigma_spn.enabled = true
					matte_Sigma_chk.enabled = true
				
					if matte_Sigma_Tex != undefined then (
						matte_Sigma_Tex.enabled = true
						
					)
				)
			
			if alphaImage != undefined do (				
					-- try to update the buttons text
					try(alphaMap_btn.text = (filenameFromPath alphaImage.PBRT_BitMap.filename)
						alphaMap_btn.tooltip = alphaImage.PBRT_BitMap.filename)catch()
			)
			
			if sigmaMap != undefined do (				
					-- try to update the buttons text
					try(matte_Sigma_Tex.text = (filenameFromPath sigmaMap.PBRT_BitMap.filename)
						matte_Sigma_Tex.tooltip = sigmaMap.PBRT_BitMap.filename)catch()
			)
			
			if bumpMap != undefined do (				
					-- try to update the buttons text
					try(matte_Bumpmap.text = (filenameFromPath bumpMap.PBRT_BitMap.filename)
						matte_Bumpmap.tooltip = bumpMap.PBRT_BitMap.filename)catch()
			)
			
			if LightSamplesEnabled_chk.checked == true then (
				LightColor_col.enabled = true
--				LightSamplesTexture_btn.enabled = true
				LightMult_spn.enabled = true
				LightSamples_spn.enabled = true
			)
		)
		
		on matte_Tex_chk changed state do (
			if matte_Tex_chk.checked then matte_Texture_btn.enabled = true else matte_Texture_btn.enabled = false
		)

		on matte_Sigma_chk changed state do (
			if matte_Sigma_chk.checked then matte_Sigma_Tex.enabled = true else matte_Sigma_Tex.enabled = false
		)

		on matte_Bump_chk changed state do (
			if matte_Bump_chk.checked then matte_Bumpmap.enabled = true else matte_Bumpmap.enabled = false
		)

		on LightSamplesEnabled_chk changed state do (
				LightColor_col.enabled = state
--				LightSamplesTexture_btn.enabled = state
				LightMult_spn.enabled = state
				LightSamples_spn.enabled = state
			)
		
		on Preview pressed do PBRT_MakePreview this Preview_Type.state Preview_Quality.state
			
		on Preview_Type changed state do RendMode = state
		on Preview_Quality changed state do PrevQuality = state
			
			
		on matte_GetMapColors pressed do (
				Local matList = PBRT_getMaterialList meditMaterials[activeMeditSlot]

				for m = 1 to matList.count do 
					(	
						if matList[m].KdTexture != undefined then matList[m].Kd = PBRT_FindAveColor matList[m].KdTexture	
					)
			)

		on matte_Mode changed state do (
			if matte_Mode.state == 1 then ( -- Lambertian
				matte_Sigma_spn.enabled = false
				matte_Sigma_Tex.enabled = false
				SigmaTmp = matte_Sigma_spn.value
				matte_Sigma_spn.value = 0.0
				)else (	-- Oren-Nayar
				matte_Sigma_spn.enabled = true
				matte_Sigma_Tex.enabled = true	
				if SigmaTmp != 0 then matte_Sigma_spn.value = SigmaTmp
				)
		)
		
		on matte_Texture_btn changed state do (
			
		)
		
		on params reload do (

		)
		

	)
	

	
	on create  do (
		)

	on load do (

		)
	
	on clone old do (
		this.name = uniquename this.name
		delegate.name = this.name 
	)
)

plugin material PBRT_Plastic name:"PBRT Plastic Mtl" classID:#(0x68dba118, 0x585962d9) extends:Standard replaceUI:true version:1
( 
	
	parameters main rollout:params
	(
		Kd type:#color default:[140,140,140] ui:plastic_Kd_col 
		KdTexture type:#textureMap ui:plastic_Kd_btn
		KdTextureEnable type:#boolean default:false ui:plastic_Kd_chk
		
		alphaAmount type:#float default:1.0 ui:alpha_spn
		alphaEnable type:#boolean default:false ui:alphaEnabled_chk
		alphaImage type:#textureMap ui:alphaMap_btn
		
		Ks type:#color default:[0,0,0] ui:plastic_Ks_col  
		KsTexture type:#textureMap ui:plastic_Ks_btn
		KsTextureEnable type:#boolean default:false ui:plastic_Ks_chk		
		
		Roughness type:#float default:0.0 ui:plastic_Roughness_spn
		RoughnessTexture type:#textureMap ui:plastic_Roughness_btn
		RoughnessTextureEnable type:#boolean default:false ui:plastic_Roughness_chk		
		
		bumpMap type:#textureMap ui:plastic_Bump_btn
		bumpMapEnable type:#boolean default:false ui:plastic_Bump_chk
		bumpMapAmount type:#float default:0.05 ui:plastic_Bump_spn
		
		on Kd set val do delegate.diffuse = val
		on KdTexture set val do delegate.diffuseMap = val	
		on KdTextureEnable set val do delegate.diffuseMapEnable = val
				
		on alphaAmount set val do delegate.opacityMapAmount = val	
		on alphaEnable set val do delegate.opacityMapEnable = val	
		on alphaImage set val do delegate.opacityMap = val		
			
		on Ks set val do delegate.specularLevel = val.r
		on KsTexture set val do delegate.specularMap = val	
		on KsTextureEnable set val do delegate.specularLevelMapEnable = val
			
		on Roughness set val do delegate.glossiness = val
		on RoughnessTexture set val do delegate.specularLevelMap = val
		on RoughnessTextureEnable set val do delegate.specularLevelMapEnable = val
			
		on bumpMap set val do delegate.bumpMap = val
		on bumpMapAmount set val do delegate.bumpMapAmount = val
		on bumpMapEnable set val do delegate.bumpMapEnable = val
			
	)

	
	rollout params "Plastic Parameters"
	(
		GroupBox Group_plastic_Kd "Kd Diffuse Reflectivity" pos:[5,0] width:315 height:94 
		colorPicker plastic_Kd_col "Diffuse Color: " pos:[12,15] width:139 height:20
		button plastic_GetMapColors "Update Color from Bitmap" pos:[171,17] width:130 height:17 enabled:false toolTip:"Get Diffuse Color from Bitmap"
		label plastic_Kd_lbl "Diffuse Texture:" pos:[12,44] width:80 height:20
		checkbox plastic_Kd_chk "" pos:[98,44] width:15 height:15 checked:false	
		mapButton plastic_Kd_btn "None" pos:[120,42] width:195 height:20 
		
		spinner alpha_spn "Alpha:" pos:[32,71] width:58 height:16 range:[0,1,1] type:#integer
		checkbox alphaEnabled_chk "" pos:[98,71] width:24 height:15
		mapButton alphaMap_btn "Alpah Map" pos:[118,69] width:195 height:20		

		GroupBox Group_plastic_Ks "Ks Specular Reflectivity" pos:[5,96] width:315 height:43	
		colorPicker plastic_Ks_col "" pos:[13,115] width:60 height:20
		checkbox plastic_Ks_chk "" pos:[89,115] width:15 height:15 checked:false
		mapButton plastic_Ks_btn "None" pos:[117,112] width:198 height:20 						
		
		GroupBox Group_Roughness "Roughness" pos:[5,143] width:315 height:44				
		spinner plastic_Roughness_spn "" pos:[13,163] width:70 height:16 range:[0,1,0.1]
		checkbox plastic_Roughness_chk "" pos:[95,165] width:15 height:15 checked:false
		mapButton plastic_Roughness_btn "None" pos:[117,161] width:198 height:20 		
		
		GroupBox plastic_Bump "Bump" pos:[5,188] width:315 height:40		
		spinner plastic_Bump_spn "" pos:[13,206] width:70 height:16 range:[-1,1,0.05]
		checkbox plastic_Bump_chk "" pos:[94,206] width:15 height:15 checked:true	
		mapButton plastic_Bump_btn "None" pos:[117,203] width:198 height:20
	
		GroupBox grp_Preview "" pos:[6,224] width:315 height:75						
		radiobuttons Preview_Obj "" pos:[15,237] width:209 height:16 labels:#("Box", "Sphere", "Torus", "Plane", "Surface") default:2 columns:5
		radiobuttons Preview_Type "" pos:[100,260] width:209 height:16 labels:#("Direct Lighting", "Metropolis") default:1 columns:2 
		radiobuttons Preview_Quality "" pos:[100,280] width:153 height:16 labels:#("Low", "Med", "High", "Best") default:2 columns:4
		button preview "Preview" pos:[15,260] width:62 height:35
		
		on params open do (
			if KdTexture != undefined then plastic_GetMapColors.enabled = true
			
			-- try to update the buttons text
			try(plastic_Kd_btn.text = (filenameFromPath KdTexture.PBRT_BitMap.filename)
				plastic_Kd_btn.tooltip = KdTexture.PBRT_BitMap.filename)catch()
				
			try(alphaMap_btn.text = (filenameFromPath alphaImage.PBRT_BitMap.filename)
				alphaMap_btn.tooltip = alphaImage.PBRT_BitMap.filename)catch()	
								
			try(plastic_Ks_btn.text = (filenameFromPath KsTexture.PBRT_BitMap.filename)
				plastic_Ks_btn.tooltip = KsTexture.PBRT_BitMap.filename)catch()
								
			try(plastic_Roughness_btn.text = (filenameFromPath KsTexture.PBRT_BitMap.filename)
				plastic_Roughness_btn.tooltip = KsTexture.PBRT_BitMap.filename)catch()	
												
			try(plastic_Bump_btn.text = (filenameFromPath bumpMap.PBRT_BitMap.filename)
				plastic_Bump_btn.tooltip = bumpMap.PBRT_BitMap.filename)catch()	
		)
		
		on Preview pressed do PBRT_MakePreview this Preview_Type.state Preview_Quality.state
	
			
		on plastic_Kd_chk changed state do (
		--	if plastic_Kd_chk.checked then plastic_Kd_btn.enabled = true else matte_Texture_btn.enabled = false
		)

		on plastic_chk changed state do (
		--	if matte_Bump_chk.checked then matte_Bumpmap.enabled = true else matte_Bumpmap.enabled = false
		)

		
		on plastic_GetMapColors pressed do (
				Local matList = PBRT_getMaterialList meditMaterials[activeMeditSlot]
	
				for m = 1 to matList.count do 
					(	
						if matList[m].KdTexture != undefined then matList[m].Kd = PBRT_FindAveColor matList[m].KdTexture	
					)
			)
			
	)
	
	
	
	on create  do (
		)

	on load do (
	--	try(if KdTexture != undefined then plastic_GetMapColors.enabled = true)catch()
	)
	
	on clone old do (
		this.name = uniquename this.name
		delegate.name = this.name 

	)
)

plugin material PBRT_Translucent name:"PBRT Translucent Mtl" classID:#(0x6a31b5d3, 0x682fa674) extends:Standard replaceUI:true version:1
( 
	
	parameters main rollout:params
	(
		Kd type:#color default:[140,140,140] ui:Translucent_Kd_col 
		KdTexture type:#textureMap ui:Translucent_Kd_btn
		KdTextureEnable type:#boolean default:false ui:Translucent_Kd_chk
		
		alphaAmount type:#float default:1.0 ui:alpha_spn
		alphaEnable type:#boolean default:false ui:alphaEnabled_chk
		alphaImage type:#textureMap ui:alphaMap_btn
		
		Ks type:#color default:[0,0,0] ui:Translucent_Ks_col  
		KsTexture type:#textureMap ui:Translucent_Ks_btn
		KsTextureEnable type:#boolean default:false ui:Translucent_Ks_chk		

		reflect type:#float default:0.5 ui:Translucent_reflect_spn  
		reflectTexture type:#textureMap ui:Translucent_reflect_btn
		reflectTextureEnable type:#boolean default:false ui:Translucent_reflect_chk	

		transmit type:#float default:0.5 ui:Translucent_transmit_spn  
		transmitTexture type:#textureMap ui:Translucent_transmit_btn
		transmitTextureEnable type:#boolean default:false ui:Translucent_transmit_chk	
		
		Roughness type:#float default:0.0 ui:Translucent_Roughness_spn
		RoughnessTexture type:#textureMap ui:Translucent_Roughness_btn
		RoughnessTextureEnable type:#boolean default:false ui:Translucent_Roughness_chk		
		
		bumpMap type:#textureMap ui:Translucent_Bump_btn
		bumpMapEnable type:#boolean default:false ui:Translucent_Bump_chk
		bumpMapAmount type:#float default:0.05 ui:Translucent_Bump_spn
		
		on Kd set val do delegate.diffuse = val
		on KdTexture set val do delegate.diffuseMap = val	
		on KdTextureEnable set val do delegate.diffuseMapEnable = val
	
		on Ks set val do delegate.specularLevel = val.r
		on KsTexture set val do delegate.specularMap = val	
		on KsTextureEnable set val do delegate.specularLevelMapEnable = val
			
		on Roughness set val do delegate.glossiness = val
		on RoughnessTexture set val do delegate.specularLevelMap = val
		on RoughnessTextureEnable set val do delegate.specularLevelMapEnable = val
			
		on bumpMap set val do delegate.bumpMap = val
		on bumpMapAmount set val do delegate.bumpMapAmount = val
		on bumpMapEnable set val do delegate.bumpMapEnable = val
			
	)

	
	rollout params "Matte Parameters"
	(
		GroupBox Group_Translucent_Kd "Kd Diffuse Reflectivity" pos:[5,0] width:315 height:94 
		colorPicker Translucent_Kd_col "Diffuse Color: " pos:[12,15] width:139 height:20
		button Translucent_GetMapColors "Update Color from Bitmap" pos:[171,17] width:130 height:17 enabled:false toolTip:"Get Diffuse Color from Bitmap"
		label Translucent_Kd_lbl "Diffuse Texture:" pos:[12,44] width:80 height:20
		checkbox Translucent_Kd_chk "" pos:[98,44] width:15 height:15 checked:false	
		mapButton Translucent_Kd_btn "None" pos:[120,42] width:195 height:20 
		
		spinner alpha_spn "Alpha:" pos:[32,71] width:58 height:16 range:[0,1,1] type:#integer
		checkbox alphaEnabled_chk "" pos:[98,71] width:24 height:15
		mapButton alphaMap_btn "None" pos:[119,69] width:195 height:20	
		
		GroupBox Group_Translucent_Ks "Ks Specular Reflectivity" pos:[5,96] width:315 height:43	
		colorPicker Translucent_Ks_col "" pos:[13,113] width:60 height:20
		checkbox Translucent_Ks_chk "" pos:[89,115] width:15 height:15 checked:false
		mapButton Translucent_Ks_btn "None" pos:[117,112] width:198 height:20 												

		GroupBox Group_Reflect "Fraction of Light Reflected" pos:[6,140] width:315 height:43
		checkbox Translucent_reflect_chk "" pos:[90,159] width:15 height:15 checked:false
		mapButton Translucent_reflect_btn "None" pos:[118,156] width:198 height:20
		spinner Translucent_reflect_spn "" pos:[11,160] width:70 height:16 range:[0,1,0.1]
		
		GroupBox grp9 "Fraction of Light Transmited" pos:[6,185] width:315 height:43
		checkbox Translucent_transmit_chk "" pos:[90,204] width:15 height:15 checked:false
		mapButton Translucent_transmit_btn "None" pos:[118,201] width:198 height:20
		spinner Translucent_transmit_spn "" pos:[11,205] width:70 height:16 range:[0,1,0.1]
		
		GroupBox Group_Roughness "Roughness" pos:[6,230] width:315 height:44				
		spinner Translucent_Roughness_spn "" pos:[13,250] width:70 height:16 range:[0,1,0.1]
		checkbox Translucent_Roughness_chk "" pos:[96,252] width:15 height:15 checked:false
		mapButton Translucent_Roughness_btn "None" pos:[118,248] width:198 height:20 		
		
		GroupBox Translucent_Bump "Bump" pos:[7,275] width:315 height:40		
		spinner Translucent_Bump_spn "" pos:[15,293] width:70 height:16 range:[-1,1,0.05]
		checkbox Translucent_Bump_chk "" pos:[96,293] width:15 height:15 checked:true		
		mapButton Translucent_Bump_btn "None" pos:[119,290] width:198 height:20
		
		GroupBox grp_Preview "" pos:[6,314] width:315 height:75				
		
		radiobuttons Preview_Obj "" pos:[15,327] width:209 height:16 labels:#("Box", "Sphere", "Torus", "Plane", "Surface") default:2 columns:5
		radiobuttons Preview_Type "" pos:[100,350] width:209 height:16 labels:#("Direct Lighting", "Metropolis") default:1 columns:2 
		radiobuttons Preview_Quality "" pos:[100,370] width:153 height:16 labels:#("Low", "Med", "High", "Best") default:2 columns:4
		button preview "Preview" pos:[15,350] width:62 height:35
		
		on params open do (				
	
			-- try to update the buttons text				

			if KdTexture != undefined do (	
					Translucent_GetMapColors.enabled = true
				
					try(Translucent_Kd_btn.text = (filenameFromPath KdTexture.PBRT_BitMap.filename)
						Translucent_Kd_btn.tooltip = KdTexture.PBRT_BitMap.filename)catch()
					) 

			if alphaImage != undefined do (				
					try(alphaMap_btn.text = (filenameFromPath alphaImage.PBRT_BitMap.filename)
						alphaMap_btn.tooltip = alphaImage.PBRT_BitMap.filename)catch()
					) 

			if KsTexture != undefined do (				
					try(Translucent_Ks_btn.text = (filenameFromPath KsTexture.PBRT_BitMap.filename)
						Translucent_Ks_btn.tooltip = KsTexture.PBRT_BitMap.filename)catch()
					) 

			if reflectTexture != undefined do (				
					try(Translucent_reflect_btn.text = (filenameFromPath reflectTexture.PBRT_BitMap.filename)
						Translucent_reflect_btn.tooltip = reflectTexture.PBRT_BitMap.filename)catch()
					) 

			if transmitTexture != undefined do (				
					try(Translucent_transmit_btn.text = (filenameFromPath transmitTexture.PBRT_BitMap.filename)
						Translucent_transmit_btn.tooltip = transmitTexture.PBRT_BitMap.filename)catch()
					)

			if RoughnessTexture != undefined do (				
					try(Translucent_Roughness_btn.text = (filenameFromPath RoughnessTexture.PBRT_BitMap.filename)
						Translucent_Roughness_btn.tooltip = RoughnessTexture.PBRT_BitMap.filename)catch()
					) 	
		
			if bumpMap != undefined do (				
					try(Translucent_Bumpmap.text = (filenameFromPath bumpMap.PBRT_BitMap.filename)
						Translucent_Bumpmap.tooltip = bumpMap.PBRT_BitMap.filename)catch()
					)
		)
		
		on Preview pressed do PBRT_MakePreview this Preview_Type.state Preview_Quality.state

		
		on Translucent_GetMapColors pressed do (
				getMaterialList meditMaterials[activeMeditSlot]
	
				for m = 1 to matList.count do 
				(
				SetDiffuseColor matList[m]	
				)
			)
			
	)
	
	
	
	on create  do (
		)

	on load do (
		try(if KdTexture != undefined then Translucent_GetMapColors.enabled = true)catch()
	)
	
	on clone old do (
		this.name = uniquename this.name
		delegate.name = this.name 

	)
)

plugin material PBRT_Substrate name:"PBRT Substrate Mtl" classID:#(0x3bf3aced, 0x71fc45bf) extends:Standard replaceUI:true version:1
( 
	
	parameters main rollout:params
	(
		Kd type:#color default:[140,140,140] ui:Substrate_Kd_col 
		KdTexture type:#textureMap ui:Substrate_Kd_btn
		KdTextureEnable type:#boolean default:false ui:Substrate_Kd_chk
		
		Ks type:#color default:[0,0,0] ui:Substrate_Ks_col  
		KsTexture type:#textureMap ui:Substrate_Ks_btn
		KsTextureEnable type:#boolean default:false ui:Substrate_Ks_chk		
		
		uRoughness type:#float default:0.0 ui:Substrate_U_Roughness_spn
		uRoughnessTexture type:#textureMap ui:Substrate_U_Roughness_btn
		uRoughnessTextureEnable type:#boolean default:false ui:Substrate_U_Roughness_chk		
		
		vRoughness type:#float default:0.0 ui:Substrate_V_Roughness_spn
		vRoughnessTexture type:#textureMap ui:Substrate_V_Roughness_btn
		vRoughnessTextureEnable type:#boolean default:false ui:Substrate_V_Roughness_chk	
		
		bumpMap type:#textureMap ui:Substrate_Bump_btn
		bumpMapEnable type:#boolean default:false ui:Substrate_Bump_chk
		bumpMapAmount type:#float default:0.05 ui:Substrate_Bump_spn
		
		on Kd set val do delegate.diffuse = val
		on KdTexture set val do delegate.diffuseMap = val	
		on KdTextureEnable set val do delegate.diffuseMapEnable = val
	
		on Ks set val do delegate.specularLevel = val.r
		on KsTexture set val do delegate.specularMap = val	
		on KsTextureEnable set val do delegate.specularLevelMapEnable = val
			
		on bumpMap set val do delegate.bumpMap = val
		on bumpMapAmount set val do delegate.bumpMapAmount = val
		on bumpMapEnable set val do delegate.bumpMapEnable = val
			
	)

	
	rollout params "Substrate Parameters"
	(
		GroupBox Group_plastic_Kd "Kd Diffuse Reflectivity" pos:[5,0] width:315 height:66 
		colorPicker Substrate_Kd_col "Diffuse Color: " pos:[12,15] width:139 height:20
		button Substrate_GetMapColors "Update Color from Bitmap" pos:[171,17] width:130 height:17 enabled:false toolTip:"Get Diffuse Color from Bitmap"
		label Substrate_Kd_lbl "Diffuse Texture:" pos:[12,44] width:80 height:20
		checkbox Substrate_Kd_chk "" pos:[98,44] width:15 height:15 checked:false	
		mapButton Substrate_Kd_btn "None" pos:[120,42] width:195 height:20 

		GroupBox Group_plastic_Ks "Ks Specular Reflectivity" pos:[5,70] width:315 height:43	
		colorPicker Substrate_Ks_col "" pos:[13,87] width:60 height:20
		checkbox Substrate_Ks_chk "" pos:[89,89] width:15 height:15 checked:false
		mapButton Substrate_Ks_btn "None" pos:[117,86] width:198 height:20 						
		
		GroupBox Group_Roughness "Roughness" pos:[5,117] width:315 height:67	
		
		spinner Substrate_U_Roughness_spn "U:" pos:[20,137] width:60 height:16 range:[0,1,0.1]
		checkbox Substrate_U_Roughness_chk "" pos:[95,139] width:15 height:15 checked:false
		mapButton Substrate_U_Roughness_btn "None" pos:[117,135] width:198 height:20 	
		spinner Substrate_V_Roughness_spn "V:" pos:[20,162] width:60 height:16 range:[0,1,0.1]
		checkbox Substrate_V_Roughness_chk "" pos:[96,164] width:15 height:15 checked:false
		mapButton Substrate_V_Roughness_btn "None" pos:[118,160] width:198 height:20	
		
		GroupBox Substrate_Bump "Bump" pos:[5,186] width:315 height:40		
		spinner Substrate_Bump_spn "" pos:[13,204] width:70 height:16 range:[-1,1,0.05]
		checkbox Substrate_Bump_chk "" pos:[94,204] width:15 height:15 checked:true		
		mapButton Substrate_Bump_btn "None" pos:[117,201] width:198 height:20

		GroupBox grp_Preview "" pos:[6,224] width:315 height:75				
		
		radiobuttons Preview_Obj "" pos:[15,237] width:209 height:16 labels:#("Box", "Sphere", "Torus", "Plane", "Surface") default:2 columns:5
		radiobuttons Preview_Type "" pos:[100,260] width:209 height:16 labels:#("Direct Lighting", "Metropolis") default:1 columns:2 
		radiobuttons Preview_Quality "" pos:[100,280] width:153 height:16 labels:#("Low", "Med", "High", "Best") default:2 columns:4
		button preview "Preview" pos:[15,260] width:62 height:35
	
	
		on params open do (
			
			if KdTexture != undefined do (	
					Substrate_GetMapColors.enabled = true
				
					try(Substrate_Kd_btn.text = (filenameFromPath KdTexture.PBRT_BitMap.filename)
						Substrate_Kd_btn.tooltip = KdTexture.PBRT_BitMap.filename)catch()
					)  
					
			if KsTexture != undefined do (				
					try(Substrate_Ks_btn.text = (filenameFromPath KsTexture.PBRT_BitMap.filename)
						Substrate_Ks_btn.tooltip = KsTexture.PBRT_BitMap.filename)catch()
					) 

			if uRoughnessTexture != undefined do (				
					try(Substrate_U_Roughness_btn.text = (filenameFromPath uRoughnessTexture.PBRT_BitMap.filename)
						Substrate_U_Roughness_btn.tooltip = uRoughnessTexture.PBRT_BitMap.filename)catch()
					) 	

			if vRoughnessTexture != undefined do (				
					try(Substrate_V_Roughness_btn.text = (filenameFromPath vRoughnessTexture.PBRT_BitMap.filename)
						Substrate_V_Roughness_btn.tooltip = vRoughnessTexture.PBRT_BitMap.filename)catch()
					) 	
					
			if bumpMap != undefined do (				
					try(Translucent_Bumpmap.text = (filenameFromPath bumpMap.PBRT_BitMap.filename)
						Translucent_Bumpmap.tooltip = bumpMap.PBRT_BitMap.filename)catch()
					)
		)
		
		on Preview pressed do PBRT_MakePreview this Preview_Type.state Preview_Quality.state
			
		
		on plastic_GetMapColors pressed do (
				getMaterialList meditMaterials[activeMeditSlot]
	
				for m = 1 to matList.count do 
				(
				SetDiffuseColor matList[m]	
				)
			)
			
	)
	
	
	
	on create  do (
		)

	on load do (
		try(if KdTexture != undefined then Substrate_GetMapColors.enabled = true)catch()
	)
	
	on clone old do (
		this.name = uniquename this.name
		delegate.name = this.name 

	)
)

plugin material PBRT_Uber name:"PBRT Uber Mtl" classID:#(0x4ea424ac, 0x5d6f8978) extends:Standard replaceUI:true version:1
( 	
	
	parameters main rollout:params
	(
		Kd type:#color default:[140,140,140] ui:uber_Kd_col 
		KdTexture type:#textureMap ui:uber_Kd_btn
		KdTextureEnable type:#boolean default:false ui:uber_Kd_chk
		
		Ks type:#color default:[0,0,0] ui:uber_Ks_col  
		KsTexture type:#textureMap ui:uber_Ks_btn
		KsTextureEnable type:#boolean default:false ui:uber_Ks_chk		
		
		Kr type:#color default:[0,0,0] ui:uber_Kr_col  
		KrTexture type:#textureMap ui:uber_Kr_btn
		KrTextureEnable type:#boolean default:false ui:uber_Kr_chk
		
		opacity type:#color default:[255,255,255] ui:uber_Opacity_col
		opacityTexture type:#textureMap ui:uber_Opacity_btn
		opacityTextureEnable type:#boolean default:false ui:uber_Opacity_chk
		
		Roughness type:#float default:0.0 ui:uber_Roughness_spn
		RoughnessTexture type:#textureMap ui:uber_Roughness_btn
		RoughnessTextureEnable type:#boolean default:false ui:uber_Roughness_chk		
		
		bumpMap type:#textureMap ui:uber_Bump_btn
		bumpMapEnable type:#boolean default:false ui:uber_Bump_chk
		bumpMapAmount type:#float default:0.05 ui:uber_Bump_spn
		
		LightColorEnable type:#boolean default:false ui:LightSamplesEnabled_chk
		emitterInv type:#boolean default:false ui:LightInvisible_chk
		LightColor type:#color default:[255,255,255] ui:lightColor_col
		LightMult type:#float default:1.0 ui:LightMult_spn
		LightSamples type:#integer default:4 ui:LightSamples_spn
		
		on Kd set val do delegate.diffuse = val
		on KdTexture set val do delegate.diffuseMap = val	
		on KdTextureEnable set val do delegate.diffuseMapEnable = val
	
		on Ks set val do delegate.specularLevel = val.r
		on KsTexture set val do delegate.specularMap = val	
		on KsTextureEnable set val do delegate.specularMapEnable = val
	
		on Kr set val do delegate.specularLevel = val.r
		on KrTexture set val do delegate.specularLevelMap = val	
		on KrTextureEnable set val do delegate.specularLevelMapEnable = val
			
		on opacity set val do (
				delegate.opacity = PBRT_clampRange stVal:0 edVal:255 nStVal:0.0 nEdVal:100.0 val:val.r
				delegate.FilterColor = val
			)
		on opacityTexture set val do delegate.opacityMap = val	
		on opacityTextureEnable set val do delegate.opacityMapEnable = val
			
		on Roughness set val do delegate.glossiness = val
		on RoughnessTexture set val do delegate.glossinessMap = val
		on RoughnessTextureEnable set val do delegate.glossinessMapEnable = val
			
		on bumpMap set val do delegate.bumpMap = val
		on bumpMapAmount set val do delegate.bumpMapAmount = val
		on bumpMapEnable set val do delegate.bumpMapEnable = val
			
		on LightColorEnable set val do delegate.useSelfIllumColor = val		
		on LightColor set val do delegate.selfIllumColor = val
		on LightMult set val do delegate.selfillumMapAmount = val			
	)

	
	rollout params "Uber Parameters"
	(
		GroupBox 		Group_uber_Kd "Kd Diffuse Reflection" pos:[5,0] width:315 height:66 
		colorPicker 	uber_Kd_col "Diffuse Color: " pos:[12,15] width:139 height:20
		button 			uber_GetMapColors "Update Color from Bitmap" pos:[171,17] width:130 height:17 enabled:false toolTip:"Get Diffuse Color from Bitmap"
		label 			uber_Kd_lbl "Diffuse Texture:" pos:[12,44] width:80 height:20
		checkbox 		uber_Kd_chk "" pos:[98,44] width:15 height:15 checked:false	
		mapButton 	uber_Kd_btn "None" pos:[120,42] width:195 height:20 

		GroupBox 		Group_uber_Ks "Ks Glossy Reflection" pos:[5,70] width:315 height:43	
		colorPicker 	uber_Ks_col "" pos:[13,87] width:60 height:20
		checkbox 		uber_Ks_chk "" pos:[89,89] width:15 height:15 checked:false
		mapButton 	uber_Ks_btn "None" pos:[117,86] width:198 height:20 						
		
		GroupBox 		Group_uber_Kr "Kr Specular Reflection" pos:[5,115] width:315 height:43	
		colorPicker 	uber_Kr_col "" pos:[13,132] width:60 height:20
		checkbox 		uber_Kr_chk "" pos:[89,134] width:15 height:15 checked:false
		mapButton 	uber_Kr_btn "None" pos:[117,131] width:198 height:20 		
		
		GroupBox 		Group_Roughness "Roughness" pos:[5,161] width:315 height:44				
		spinner 		uber_Roughness_spn "" pos:[13,181] width:70 height:16 range:[0,1,0.1]
		checkbox 		uber_Roughness_chk "" pos:[95,183] width:15 height:15 checked:false
		mapButton 	uber_Roughness_btn "None" pos:[117,179] width:198 height:20 		
		
		GroupBox 		Group_uber_Opacity "Opacity" pos:[6,207] width:315 height:43
		colorPicker 	uber_Opacity_col "" pos:[14,224] width:60 height:20 color:(color 0 0 155)
		checkbox 		uber_Opacity_chk "" pos:[90,226] width:15 height:15 checked:false
		mapButton 	uber_Opacity_btn "None" pos:[118,223] width:198 height:20	
		
		GroupBox 		uber_Bump "Bump" pos:[5,254] width:315 height:40		
		spinner 		uber_Bump_spn "" pos:[13,272] width:70 height:16 range:[-1,1,0.05]
		checkbox 		uber_Bump_chk "" pos:[94,272] width:15 height:15 checked:true	
		mapButton 	uber_Bump_btn "None" pos:[117,269] width:198 height:20	
	
	
		GroupBox LightSamples_gup "Emitter" pos:[6,294] width:315 height:67
		checkbox LightSamplesEnabled_chk "On/Off:" pos:[18,311] width:57 height:13 	
		checkbox LightInvisible_chk "Invisible:" pos:[90,311] width:57 height:13		
		colorPicker lightColor_col "Light Color" pos:[180,309] width:92 height:20 enabled:false
		spinner LightMult_spn "Light Multiplier:" pos:[70,337] width:80 height:16 enabled:false range:[0,1000,1]
		spinner LightSamples_spn "Light Samples:" pos:[210,338] width:80 height:16 enabled:false range:[0,1000,4] type:#integer	
		
		GroupBox grp_Preview "" pos:[6,364] width:315 height:75				
		
		radiobuttons Preview_Obj "" pos:[15,378] width:209 height:16 labels:#("Box", "Sphere", "Torus", "Plane", "Surface") default:2 columns:5
		radiobuttons Preview_Type "" pos:[100,400] width:209 height:16 labels:#("Direct Lighting", "Metropolis") default:1 columns:2 
		radiobuttons Preview_Quality "" pos:[100,420] width:153 height:16 labels:#("Low", "Med", "High", "Best") default:2 columns:4
		button preview "Preview" pos:[15,400] width:62 height:35


		on params open do (
			if KdTexture != undefined then uber_GetMapColors.enabled = true
			
			if LightSamplesEnabled_chk.checked == true then (
				LightColor_col.enabled = true
				LightMult_spn.enabled = true
				LightSamples_spn.enabled = true
			)			
		)	
		
		on LightSamplesEnabled_chk changed state do (
				LightColor_col.enabled = state
				LightMult_spn.enabled = state
				LightSamples_spn.enabled = state
		)
		
		on Preview pressed do PBRT_MakePreview this Preview_Type.state Preview_Quality.state

			
		on uber_Kd_chk changed state do (
		--	if plastic_Kd_chk.checked then plastic_Kd_btn.enabled = true else matte_Texture_btn.enabled = false
		)

		on uber_chk changed state do (
		--	if matte_Bump_chk.checked then matte_Bumpmap.enabled = true else matte_Bumpmap.enabled = false
		)

		
		on uber_GetMapColors pressed do (
				getMaterialList meditMaterials[activeMeditSlot]
	
				for m = 1 to matList.count do 
				(
				SetDiffuseColor matList[m]	
				)
			)
			
	)
	
	
	
	on create  do (
		)

	on load do (
		try(if KdTexture != undefined then uber_GetMapColors.enabled = true)catch()
	)
	
	on clone old do (
		this.name = uniquename this.name
		delegate.name = this.name 

	)
)

plugin material PBRT_Metal name:"PBRT Metal Mtl" classID:#(0x155342a8, 0x4b96def6) extends:Standard replaceUI:true version:1
( 
	
	parameters main rollout:params
	(
		
		eta type:#string default:"" ui:metal_eta_btn
		k type:#string ui:metal_k_btn
		
		Roughness type:#float default:0.0 ui:metal_Roughness_spn
		RoughnessTexture type:#textureMap ui:metal_Roughness_btn
		RoughnessTextureEnable type:#boolean default:false ui:metal_Roughness_chk		
		
		bumpMap type:#textureMap ui:metal_Bump_btn
		bumpMapEnable type:#boolean default:false ui:metal_Bump_chk
		bumpMapAmount type:#float default:0.05 ui:metal_Bump_spn
		
		on Roughness set val do delegate.glossiness = val
		on RoughnessTexture set val do delegate.glossinessMap = val
		on RoughnessTextureEnable set val do delegate.glossinessMapEnable = val
			
		on bumpMap set val do delegate.bumpMap = val
		on bumpMapAmount set val do delegate.bumpMapAmount = val
		on bumpMapEnable set val do delegate.bumpMapEnable = val
				
	)

	
	rollout params "Metal Parameters"
	(
		Button 			metal_eta_btn "None" pos:[40,7] width:279 height:20
		label 			metal_eta_lbl 	"Eta:" pos:[12,9] width:21 height:14
		Button 			metal_k_btn 	"None" pos:[41,33] width:279 height:20
		label 			metal_k_lbl 	"K:" pos:[13,35] width:21 height:14	
		
		GroupBox Group_Roughness "Roughness" pos:[5,57] width:315 height:44				
		spinner 	metal_Roughness_spn "" pos:[13,77] width:70 height:16 range:[0,1,0.1]
		checkbox 	metal_Roughness_chk "" pos:[95,79] width:15 height:15 checked:false
		mapButton metal_Roughness_btn "None" pos:[117,75] width:198 height:20 		
		
		GroupBox 		_metal_Bump "Bump" pos:[5,104] width:315 height:40		
		spinner 		metal_Bump_spn "" pos:[13,122] width:70 height:16 range:[-1,1,0.05]
		checkbox 		metal_Bump_chk "" pos:[94,122] width:15 height:15 checked:true	
		mapButton 	metal_Bump_btn "None" pos:[117,119] width:198 height:20	
	
		GroupBox grp_Preview "" pos:[6,144] width:315 height:75				
		
		radiobuttons Preview_Obj "" pos:[15,157] width:209 height:16 labels:#("Box", "Sphere", "Torus", "Plane", "Surface") default:2 columns:5
		radiobuttons Preview_Type "" pos:[100,180] width:209 height:16 labels:#("Direct Lighting", "Metropolis") default:1 columns:2 
		radiobuttons Preview_Quality "" pos:[100,200] width:153 height:16 labels:#("Low", "Med", "High", "Best") default:2 columns:4
		button preview "Preview" pos:[15,180] width:62 height:35
		
		on params open do (		

			-- try to update the buttons text
			try(metal_eta_btn.text = (filenameFromPath eta)
				metal_eta_btn.tooltip = eta)catch()				
			
			try(metal_k_btn.text = (filenameFromPath k)
				metal_k_btn.tooltip = k)catch()	
				
			try(metal_Bump_btn.text = (filenameFromPath bumpMap.PBRT_BitMap.filename)
					metal_Bump_btn.tooltip = bumpMap.PBRT_BitMap.filename)catch()	
		)	
		
		on metal_eta_btn pressed do (
				in_name = getOpenfileName caption:"Select Eta File:" types:"Settings File(*.spd)|*.spd"
				if in_name != undefined then (
					eta = in_name
					
					-- try to update the buttons text
					try(metal_eta_btn.text = (filenameFromPath eta)
						metal_eta_btn.tooltip = eta)catch()	
					)
			)
		
		on metal_k_btn pressed do (
				in_name = getOpenfileName caption:"Select Eta File:" types:"Settings File(*.spd)|*.spd"
				if in_name != undefined then (
					k = in_name
					-- try to update the buttons text
					try(metal_k_btn.text = (filenameFromPath k)
							metal_k_btn.tooltip = k)catch()	
						)
			)		
	
			
		on Preview pressed do PBRT_MakePreview this Preview_Type.state Preview_Quality.state
			
	)
	
	
	
	on create  do (
		this.name = uniquename "Metal_Material"
		)

	on load do (

	)
	
	on clone old do (
		this.name = uniquename this.name
		delegate.name = this.name 

	)
)

plugin material PBRT_Glass name:"PBRT Glass Mtl" classID:#(0x6d92d47d, 0x4a52b7ca) extends:Standard replaceUI:true version:1
( 
		
	parameters main rollout:params
	(
		Kr type:#color default:[255,255,255] ui:glass_Kr_col
		IOR type:#float default:1.50 ui:glass_IOR_spn		
		Kt type:#color default:[255,255,255] ui:glass_Kt_col 
		
		bumpMap type:#textureMap ui:glass_Bump_btn
		bumpMapEnable type:#boolean default:false ui:glass_Bump_chk
		bumpMapAmount type:#float default:0.05 ui:glass_Bump_spn
		
		on Kr set val do delegate.diffuse = val
		on Kt set val do  (
			delegate.opacity = PBRT_clampRange stVal:0 edVal:255 nStVal:100.0 nEdVal:0.0 val:((val.r * val.g * val. g) /3)
			delegate.FilterColor = val
		)
		on IOR set val do delegate.ior = val	
			
		on bumpMap set val do delegate.bumpMap = val
		on bumpMapAmount set val do delegate.bumpMapAmount = val
		on bumpMapEnable set val do delegate.bumpMapEnable = val
	)

	
	rollout params "Matte Parameters"
	(
		colorPicker 	glass_Kr_col 	"Reflection:" pos:[11,9] width:112 height:20
		spinner 		glass_IOR_spn "IOR" pos:[21,36] range:[0,100,1.50] width:70 height:16 type:#float
		colorPicker 	glass_Kt_col 	"Transmission:" pos:[152,10] width:112 height:20 color:(color 0 0 155)
		
		GroupBox 		glass_Bump 		"Bump" pos:[5,61] width:315 height:40			
		spinner 		glass_Bump_spn "" pos:[13,79] width:70 height:16 range:[-1,1,0.05]
		checkbox		glass_Bump_chk "" pos:[94,79] width:15 height:15 checked:true	
		mapButton 	glass_Bump_btn "None" pos:[117,76] width:198 height:20

		GroupBox grp_Preview "" pos:[6,104] width:315 height:75				
		
		radiobuttons Preview_Obj "" pos:[15,117] width:209 height:16 labels:#("Box", "Sphere", "Torus", "Plane", "Surface") default:2 columns:5
		radiobuttons Preview_Type "" pos:[100,140] width:209 height:16 labels:#("Direct Lighting", "Metropolis") default:1 columns:2 
		radiobuttons Preview_Quality "" pos:[100,160] width:153 height:16 labels:#("Low", "Med", "High", "Best") default:2 columns:4
		button preview "Preview" pos:[15,140] width:62 height:35
		
		on Preview pressed do PBRT_MakePreview this Preview_Type.state Preview_Quality.state

	)
	
	
	on create  do (
		)

	on load do (

	)
	
	on clone old do (
		this.name = uniquename this.name
		delegate.name = this.name 

	)
)

plugin material PBRT_Light name:"PBRT Light Mtl" classID:#(0x4a3aa127, 0x20efe15) extends:Standard replaceUI:true version:1
( 
	parameters main rollout:params
	(
		LightColor type:#color default:[255,255,255] ui:lightColor_col
		LightMult type:#float default:1.0 ui:LightMult_spn
		LightSamples type:#integer default:1 ui:LightSamples_spn
		
		emitterInv type:#boolean default:false ui:LightInvisible_chk
				
		on LightColor set val do (
			delegate.selfIllumColor = val
		)
	)

	rollout params "Matte Parameters"
	(
		colorpicker LightColor_col "Light Color:" pos:[10,5]
		checkbox LightInvisible_chk "Invisible:" pos:[140,8] width:57 height:13
		spinner LightMult_spn "Light Multiplier:" range:[0.1,1000000.0,1.0] pos:[9, 35] type:#float width:90 align:#left
		spinner LightSamples_spn "Light Samples:" range:[1, 1000, 8] pos:[10, 60] type:#integer width:90 align:#left

		GroupBox grp_Preview "" pos:[6,104] width:315 height:75				
		
		radiobuttons Preview_Obj "" pos:[15,117] width:209 height:16 labels:#("Box", "Sphere", "Torus", "Plane", "Surface") default:2 columns:5
		radiobuttons Preview_Type "" pos:[100,140] width:209 height:16 labels:#("Direct Lighting", "Metropolis") default:1 columns:2 
		radiobuttons Preview_Quality "" pos:[100,160] width:153 height:16 labels:#("Low", "Med", "High", "Best") default:2 columns:4
		button preview "Preview" pos:[15,140] width:62 height:35
		
		on Preview pressed do PBRT_MakePreview this Preview_Type.state Preview_Quality.state
	)
	
		on create  do (
			delegate.useSelfIllumColor = on
		)
		
)

