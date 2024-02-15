:: TODO: The --target-env=vulkan1.0 should maybe be left out (but it's working now so I don't want to touch it)
glslc -O --target-env=vulkan1.0 -ofillLDM.spv fillLDM.comp
glslc -O --target-env=vulkan1.0 -oapplyLineshapes.spv applyLineshapes.comp

pause
