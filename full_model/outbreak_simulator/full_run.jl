# Author: Cameron Zachreson
# Institution: The University of Melbourne
# Simulation code acompanying the manuscript entitled: 
# "A model-based assessment of social isolation practices for COVID-19 outbreak response in residential care facilities"
# Date released: Dec. 18, 2023


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



include("./Main_RACF_OB_v10_LD_unmitigated.jl")

main_OB()



