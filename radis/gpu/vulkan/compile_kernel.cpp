

void compile_source(const char* source, size_t* out_size){

    // Load if not initialized yet
	if (!inputLaunchConfiguration.isCompilerInitialized) {
		if (!app->configuration.isCompilerInitialized) {
			int resGlslangInitialize = glslang_initialize_process();
			if (!resGlslangInitialize) return VKFFT_ERROR_FAILED_TO_INITIALIZE;
			app->configuration.isCompilerInitialized = 1;
		}
	}

    VkFFTApplication* app = new VkFFTApplication({});
    app->configuration.loadApplicationFromString = 0;
    app->configuration.maxComputeWorkGroupCount[0];
    app->configuration.maxComputeWorkGroupCount[1];
    app->configuration.maxComputeWorkGroupCount[2];
    app->configuration.maxComputeWorkGroupSize[0];
    app->configuration.maxComputeWorkGroupSize[1];
    app->configuration.maxComputeWorkGroupSize[2];
    app->configuration.halfPrecision;

    VkFFTAxis* axis = new VkFFTAxis({});
    VkFFTSpecializationConstantsLayout specconst;
    specconst.code0 = source;
    axis->specializationConstants = specconst;

    VkFFTResult res = VK_SUCCESS;
    res = VkFFT_CompileKernel(VkFFTApplication* app, VkFFTAxis* axis)
    if (res != VKFFT_ERROR_FAILED_TO_CREATE_SHADER_MODULE){
        // Error!
    }
    else
    {
        res = VK_SUCCESS
    }

    // axis->binary = code;


    // free if initialized
    if (app->configuration.isCompilerInitialized) {
		glslang_finalize_process();
		app->configuration.isCompilerInitialized = 0;
	}
}
