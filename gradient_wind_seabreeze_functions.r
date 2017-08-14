# gradient_wind_seabreeze_functions.r

deg <- function(deg) {
	if(is.na(deg)) {
		return(NA)
	} else {
		while(deg >= 360 | deg < 0) {
			if(deg < 360) {
				deg <- deg + 360
			} else if(deg >= 360){
				deg <- deg - 360
			}
		}
		return(deg)
	}
}
