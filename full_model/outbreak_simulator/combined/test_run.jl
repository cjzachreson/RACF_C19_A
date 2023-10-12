#include("./Main_RACF_OB_v9.jl")

include("./header_RACF.jl")

include("./Setup_RACF_v9.jl")
import .Setup_RACF

include("./Networks_RACF_v9.jl")
import .Networks_RACF

include("./Diseases_RACF_v9.jl")
import .Diseases_RACF

include("./Agents_RACF_v9.jl")
import .Agents_RACF

include("./Facility_Structure_v9.jl")
import .Facility_Structure

include("./Outbreak_Response_RACF_v9.jl")
import .Outbreak_Response

include("./Transmission_Dynamics_v9.jl")
import .Transmission_Dynamics


#include("./Main_RACF_R0_v9.jl")

include("./Main_RACF_R0_Hom_v9.jl")

#include("./Main_RACF_FS_v9.jl")

main_R0_homogeneous()

main_R0()

main_Final_Size()



