glslc -O -otest_shader.spv --target-env=vulkan1.3 test_shader.comp
glslc -O -otest_shader2.spv --target-env=vulkan1.3 test_shader2.comp
glslc -O -otest_spec.spv --target-env=vulkan1.3 test_spec.comp

pause
