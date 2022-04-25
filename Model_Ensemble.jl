using Mimi
using MimiDICE2013
using MimiSNEASY
using MimiBRICK

# run DICE
DICE = MimiDICE2013.get_model()
run(DICE)

co2_emissions = getdataframe(DICE, :emissions=>:E) # E is Total CO2 emissions (GtCO2 per year)
print(co2_emissions) # prints total co2 emissions for all years of DICE

# run SNEASY-BRICK
SNEASY_BRICK = MimiBRICK.create_sneasy_brick()
run(SNEASY_BRICK)
#explore(SNEASY_BRICK)

# update parameters to feed DICE emissions to SNEASY-BRICK
#update_param!(SNEASY_BRICK, :ccm, :CO2_emissions, ??)
print(getdataframe(SNEASY_BRICK, :ccm=>:CO2_emissions))
