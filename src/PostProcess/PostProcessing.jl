# 几何体电流计算相关函数
include("CurrentOnGeos.jl")

# 辐射积分计算相关函数
include("RadiationIntegral.jl")

# 远场计算函数
include("FarField.jl")

# 场提取 (Surface Currents, etc.)
include("FieldExtraction.jl")

# 雷达散射截面(Radar Cross Section, RCS)计算相关函数
include("RCS.jl")

# S参数计算相关函数 (modular implementation)
include("SParameter/Types.jl")
include("SParameter/Conversion.jl")
include("SParameter/Excitation.jl")
include("SParameter/PortMeasurement.jl")
include("SParameter/Core.jl")
include("SParameter/Touchstone.jl")
