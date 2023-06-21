use "S:\projects\ELSAH\Analytic Tools\Test data\Gen2\TestData_May10_2023.dta", clear

cd "S:\projects\ELSAH\Analytic Tools\Stata code\reportExample"

*Keeping target ensemble + ID variables
keep *1 *ID Baseweight wave WaveResponseStatus
drop f* e*

*Renaming/Reconfiguring variables for example (excluding one-offs)
ren PreschoolID ProviderID
ren wave WaveID
ren WaveResponseStatus WaveResponse
ren Baseweight BaseWeight

ren a1 SexAtBirth 							//stable binary
ren b1 Discipline							//stable ordered categorical
ren c1 ProviderType 						//stable unordered categorical
ren d1 BaselineHouseholdIncome			 	//stable interval-valued

ren ad1 SNAPEligible 						//repeating binary
ren bd1 MotherEducation						//repeating ordered categorical
ren cd1 ProviderDelivery 					//repeating unordered categorical
recode ProviderDelivery 4 = 2 				//recode extra value to "Remote"
ren dd1 AssessmentScore 					//repeating interval-valued

*Define value labels
lab def sex ///
	0 Male  ///
	1 Female ///
	, replace

lab def ptype ///
	1 "Public school" ///
	2 "Head Start" ///
	3 "Another licensed Provider" ///
	4 "Provider not listed" ///
	, replace

lab def yesno ///
	0 "No" ///
	1 "Yes" ///
	, replace
	
lab def educ ///
	1 "Did not complete high school" ///
	2 "High school degree or equivalent" ///
	3 "Bachelor's degree" ///
	4 "Master's degree or higher" ///
	, replace

lab def pdeliver ///
	1 "In-person" ///
	2 "Remote" ///
	3 "Hybrid" ///
	, replace

*Label values
lab val SexAtBirth sex
lab val ProviderType ptype
lab val SNAPEligible yesno
lab val MotherEducation educ
lab val ProviderDelivery pdeliver

*add our Providers to stata's libraries
adopath + "S:\projects\ELSAH\Analytic Tools\Stata code\PanelPatch"

*Generating "age" variable
set seed 101 

cap drop Age
bys ChildID: gen Age = int(runiform(3,6)) if WaveID == 1
bys ChildID (WaveID): replace Age = Age[_n-1] + 1 if missing(Age)

*Creating item and wave nonreponse instances
foreach VAR in SexAtBirth Discipline ProviderType BaselineHouseholdIncome SNAPEligible MotherEducation ProviderDelivery AssessmentScore{
	
	*Creating item nonreponse
	cap drop rand
	gen rand = abs(rnormal(0,1))
	replace `VAR' = . if rand >2 //missing for ~5% of total obs
	drop rand
	
	*Enforcing wave nonresponse
	replace `VAR' = . if WaveResponse == 0
	
} //VAR

*Replacing any missing wave 1 Stable vars with nonmissing value
foreach sVAR in SexAtBirth Discipline ProviderType BaselineHouseholdIncome{
	tempvar tmpS
	qui gegen `tmpS' = max(`sVAR'), by(ChildID) //grab any nonmissing value
	bys ChildID: replace `sVAR' = `tmpS' if WaveID == 1 & missing(`sVAR') //replace wave 1 stable var
	qui drop `tmpS'
} //sVAR

*Creating Replicate weights
survwgt create jk1,psu(ProviderID) weight(BaseWeight) stem(repwgt_)

save "panelPATCHexamplePRE", replace

**
* We should provide the data set from this point.
**

*Summarize item nonresponse
missings report /// stable vars
	SexAtBirth Discipline ProviderType BaselineHouseholdIncome ///
	if WaveResponse == 1 & WaveID == 1, percent
	
missings report /// repeating vars and age
	SNAPEligible MotherEducation ProviderDelivery AssessmentScore ///
	Age ///
	if WaveResponse == 1, percent 

*Summarize wave nonreponse
tab WaveResponse

***begin imputation
PanelPatch  ///
    SexAtBirth Discipline ProviderType BaselineHouseholdIncome ///
	SNAPEligible MotherEducation ProviderDelivery AssessmentScore ///
	Age ///
    , ///
    add(5) ///
	weightvar(BaseWeight) ///
    j(ProviderID) i(ChildID) wave(WaveID) ///
    snumericvars(SexAtBirth BaselineHouseholdIncome) ///
    sunorderedvars(ProviderType) ///
    sorderedvars(Discipline)  ///
    vnumericvars(SNAPEligible AssessmentScore Age) ///
    vunorderedvars(ProviderDelivery) ///
    vorderedvars(MotherEducation) ///
    waveresponseflag(WaveResponse)
	
save "panelPATCHexamplePOST", replace