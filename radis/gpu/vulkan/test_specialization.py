import numpy as np
from vulkan_compute_lib import ArrayBuffer, ComputeApplication


def initialize():
    global app, data_in_d, data_FT_d, data_out_d, init_h
    print("Init... ", end="")

    app = ComputeApplication(deviceID=0)
    app.data_d = ArrayBuffer((10,), np.int32, binding=2, app=app)

    app.schedule_shader("test_spec.spv", (100, 1, 1), (10, 10, 1, 30))
    app.endCommandBuffer()

    print("Done!")
    return app


app = initialize()
app.run()

res = app.data_d.getData()

print(res)
