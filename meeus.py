#There are functions taken from Jean Meeus's "Astronomical Algorithms"
#There are also some helper original helper functions that did not come from Meeus, but are necessary to get a complete picture.

import math

#################################
#!!!!!REALLY IMPORTANT NOTE!!!!!##############################################
#All the angles that get passed into functions MUST be in radians            #
#I am not joking with you.                                                   #
#I wrote a conversion function called angle_convert so you can change formats#
##############################################################################

########################
#Astronomical Constants#
########################

#From IERS Technical note 32 (2004)
ae = 6.3781366e6	#Equitorial radius of the Earth (meters)
f = 1/298.26524		#Earth flattening constant (unit-less)

#Angular velocity of Earth with respect to the fixed stars, not the vernal equinox (radians/second) 1996.5 value (p83)
omega_earth = 7.292114992e-5

#Derived Constants
be = ae*(1-f)	#Polar radius of Earth (meters) (p 81-82)
def ee():	#The eccentricity of the Earth's meridian ellipse (unit-less) (p82)
 return math.sqrt(2*f-f*f)
#phi - phi_prime reaches its maximum value at u=45 degrees. Designate these phi_naught (geographical latitude) and phi_naught_prime (geocentric latitude)
phi_naught = math.atan(ae/be)
phi_naught_prime = math.atan(be/ae)

##########################
#The functions from Meeus#
##########################

#This function returns the Julian Day number since -4716 BCE (page 60)
#It takes its arguments in the Gregorian calendar
#The day may contain a decimal (see day_fraction below)

def julian_day(year, month, day):
 if month == 1 or month == 2:
  year = year - 1
  month = month + 12
 a = int(year/100)
 b = 2-a+int(a/4)
 jd = int(365.25*(year+4716))+int(30.6001*(month+1)) + b + day - 1524.5
 return jd

#The Julian Date of Jan 0.0 for a given year (page 62)
def julian_day_zero(year):
 year = year - 1
 a = int(year/100)
 jd0 = int(365.25*year)-a+int(a/4)+1721424.5
 return jd0

#The Modified Julian Day (MJD) begins at Greenwich mindnight instead of Greenwich noon.
#Often used for artifical Earth satellites (page 63)
def modified_julian_day(julian_day):
 return julian_day - 2400000.5

#Return the calendar data given the Julian Day (page 63)
#Return None for negative Julian day numbers
def calendar_date(julian_day):
 if julian_day <= 0:
  return None
 z = int(julian_day)
 f = julian_day - z
 if z < 2299161:
  a=z
 else:
  alpha = int( (z - 1867216.25)/36524.25 )
  a = z + 1 + alpha-int(alpha/4)
 b = a + 1524
 c = int( (b-122.1)/365.25 )
 d = int(365.25*c)
 e = int( (b-d)/30.6001 )
 day = b-d-int(30.6001*e)+f
 if e < 14:
  month = e - 1
 elif e == 14 or e == 15:
  month = e - 13
 if month > 2:
  year = c - 4716
 elif month == 1 or month == 2:
  year = c-4715
 return [year,month,day]

#Number of days between two calendar days (page 64)
def calendar_date_difference(year1, month1, day1, year2, month2, day2):
 jd1 = julian_day(years2,month2,day2) - julian_day(year1,month1,day1)
 return jd1

#Given the Julian Day, returns [week_day, name_of_day] (page 64)
def day_of_the_week(julian_day):
 days_of_the_week = ['Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday']
 a = (julian_day+1.5)%7
 week_day = int(a)
 return [week_day, days_of_the_week[week_day]]

#The day number for the year.  I have seen occasionally seen this referred to as the Julian Day, but that usage is incorrect (page 66)
def day_of_the_year(year,month,day):
 if is_leap_year(year):
  k=1
 else:
  k=2
 n = int( (275*month)/9) - k * int( (month+9)/12 ) + day - 30
 return n

#I think there is a typo in Meeus.  You have to add 30 to the day_number (page 66)
def date_from_day_of_the_year(year, day_number):
 if is_leap_year(year):
  k=1
 else:
  k=2
 month = int( (9*(k+day_number+30))/275 + 0.98 )
 if day_number < 32:
  month = 1
 day = day_number = int((275*month)/9) + k*int((month+9)/12)+30
 return [month,day]

#Needed for further calculations (pages 81-82)
#phi is the geographical latitude
#phi_prime is the geocentric latitude
#H is the observer's height above sea level
def rho_sin_phi_prime(phi, H):
 u = math.atan((be/ae)*math.tan(phi))
 return (be/ae)*math.tan(phi)

def rho_cos_phi_prime(phi, H):
 u = math.atan((be/ae)*math.tan(phi))
 return math.cos(u)+(H/ae)*math.cos(phi)

#Rho is the observer's distance to the center of the Earth
def rho_sea_level(phi):
 return 0.9983271+0.0016764*math.cos(2*phi)-0.0000035*math.cos(4*phi)

#The parallel of latitude phi is a circle with radius Rp. (page 83)
def Rp(phi):
 return (ae*math.cos(phi))/(math.sqrt(1-ee()*ee()*math.sin(phi)*math.sin(phi)))

#The length of one degree of longitude for a given latitude (page 83)
def longitude_degree_length(phi):
 return math.pi*Rp(phi)/180.0

#The linear velocity of a point at latitude phi due to the Earth's rotation in seconds (page 84)
def linear_velocity_at_latitude(phi):
 return omega_earth*Rp(phi)

#The radius of the curvature of the Earth's meridian at latitude phi (page 84)
def Rm(phi):
 tmp = (1-ee()*ee()*math.sin(phi)*math.sin(phi))*(1-ee()*ee()*math.sin(phi)*math.sin(phi))*(1-ee()*ee()*math.sin(phi)*math.sin(phi))
 return (ae*(1-ee()*ee()))/math.sqrt(tmp)

#One degree of latitude corresponds to a length of... (page 84)
def length_one_degree_latitude(rm):
 math.pi*rm/180.0

#Distance between two point on the geoid at sea level. Angles must be in radians (page 85)
def distance_at_sea_level(phi1,L1,phi2,L2):
 F=(phi1+phi2)/2.0
 G=(phi1-phi2)/2.0
 l=(L1-L2)/2.0
 S=math.sin(G)*math.sin(G)*math.cos(l)*math.cos(l)+math.cos(F)*math.cos(F)*math.sin(l)*math.sin(l)
 C=math.cos(G)*math.cos(G)*math.cos(l)*math.cos(l)+math.sin(F)*math.sin(F)*math.sin(l)*math.sin(l)
 l_omega = math.atan(math.sqrt(S/C))
 R=math.sqrt(S*C)/l_omega	#l_omega in radians
 D=2.0*l_omega*ae
 H1=(3*R-1)/(2*C)
 H2=(3*R+1)/(2*S)
 return D*(1+f*H1*math.sin(F)*math.sin(F)*math.cos(G)*math.cos(G)-f*H2*math.cos(F)*math.cos(F)*math.sin(G)*math.sin(G))

############
#Chapter 12#
############

#Compute the mean sidereal time at hour zero UT at greenwich (page 87)
#This is only valid if the julian_date is equal to 0h UT. All other date will give incorrect results 
def greenwich_0h_UT(julian_date):
 T=(julian_date-2451545.0)/36525
 return 100.46061837 + 36000.770053608*T + 0.000387933*T*T - T*T*T/38710000	#These values were selected by the IAU in 1982

#Compute the mean sidereal time for any instant (page 87)
def sidereal_UT(julian_date):
 theta_naught = greenwich_0h_UT(julian_date)+1.00273790935*julian_date
 return theta_naught

#Compute the mean sidereal time at Greenwich expressed in degrees (Page 87-88)
def sidereal_UT_degrees(julian_date):
 T=(julian_date-2451545.0)/36525
 return 280.46061837 + 360.98564736629*(julian_date-2451545.0) + 0.000387933*T*T - T*T*T/38710000

###############################################################################################
#BUG ALERT                                                                                    #
###############################################################################################
#Note: The apparent sidereal time requires corrections that are not explained until chapter 22#
#I will leave this corrections out until I understand them better.  This is a bug and needs to#
#be fixed.                                                                                    #
###############################################################################################

############
#Chapter 13#
############

def hour_angle_local(local_sidereal_time,right_ascension):
 return local_sidereal_time - right_ascension

def hour_angle_greenwich(greenwich_sidereal_time, observer_longitude, right_ascension):
 return greenwich_sidereal_time - observer_longitude - right_ascension

def equitorial_to_ecliptical(right_ascension, declination, obliquity_of_ecliptic):
 n=math.sin(right_ascension)*math.cos(obliquity_of_ecliptic)+math.tan(declination)*math.sin(obliquity_of_ecliptic)
 d=math.cos(right_ascension)
 ecliptical_longitude = math.atan(n/d)
 n=math.sin(declination)*math.cos(obliquity_of_ecliptic) - math.cos(declination)*math.sin(obliquity_of_ecliptic)*math.sin(right_ascension)
 ecliptical_latitude = math.asin(n)
 return [ecliptical_latitude,ecliptical_longitude]

def ecliptical_to_equatorial(ecliptical_latitude, ecliptical_longitude, obliquity_of_ecliptic):
 n=math.sin(ecliptical_longitude)*math.cos(obliquity_of_ecliptic) - math.tan(ecliptical_latitude)*math.sin(obliquity_of_ecliptic)
 d=math.cos(ecliptical_longitude)
 right_ascension=math.atan(n/d)
 n=math.sin(ecliptical_latitude)*math.cos(obliquity_of_ecliptic) + math.cos(ecliptical_latitude)*math.sin(obliquity_of_ecliptic)*math.sin(ecliptical_longitude)
 declination=math.asin(n)
 return [right_ascension, declination]

#Page 92
#Astronomers take south to be zero degrees azimuth.  This is the opposite of navigators and pretty much everyone else.
#To set the azimuth to north-as-zero-degrees, add 180 degrees to the answer.
#Add a function to do this same function with north-as-zero
def local_horizontal_coordinates_astro(hour_angle, observer_latitude, declination):
 n=math.sin(hour_angle)
 d=math.cos(hour_angle)*math.sin(observer_latitude) - math.tan(declination)*math.cos(observer_latitude)
 azimuth = math.atan(n/d)
 n=math.sin(observer_latitude)*math.sin(declination) + math.cos(observer_latitude)*math.cos(declination)*math.cos(hour_angle)
 altitude = math.asin(n)
 return [azimuth,altitude]

def local_horizontal_coordinates_nav(hour_angle, observer_latitude, declination):
 n=math.sin(hour_angle)
 d=math.cos(hour_angle)*math.sin(observer_latitude) - math.tan(declination)*math.cos(observer_latitude)
 azimuth = math.atan(n/d)
 n=math.sin(observer_latitude)*math.sin(declination) + math.cos(observer_latitude)*math.cos(declination)*math.cos(hour_angle)
 altitude = math.asin(n)
 return [azimuth+math.pi(),altitude]

#Page 94
def equitorial_to_local_horizontal(azimuth,altitude,observer_latitude):
 n=math.sin(azimuth)
 d=math.cos(azimuth)*math.sin(observer_latitude) + math.tan(altitude)*math.cos(observer_latitude)
 local_hour_angle = math.atan(n/d)
 n=math.sin(observer_latitude)*math.sin(altitude) - math.cos(observer_latitude)*math.cos(altitude)*math.cos(azimuth)
 declination = math.asin(n)
 return [local_hour_angle,declination]

#Page 94
#The galactic pole is fixed by convention for a given equinox.  Each equinox will list an right ascension and declination for the year.
###def equitorial_to_galactic(
###come back and do this one after you have covered chapter 21.  I am uncertain how to do this correctly while accounting for arbitrary epochs at this time

############
#Chapter 14#
############

#Page 98
def parallatic_angle(observer_latitude,declination_of_body,hour_angle_of_body):
 n=math.sin(hour_angle_of_body)
 d=math.tan(observer_latitude)*math.cos(declination_of_body) - math.sin(declination_of_body)*math.cos(hour_angle_of_body)
 try:
  q=math.atan(n/d)
  return q
 except:
  return 'undefined zenith'

#Page 99
#Use this when the body is on the horizon, i.e. when it is rising or setting
def parallatic_angle_horizon(observer_latitude, declination_of_body):
 n=math.sin(observer_lattitude)
 d=math.cos(declination_of_body)
 q=acos(n/d)
 return q

#Page 99
def ecliptic_longitude_at_horizon(obliquity_of_ecliptic, observer_latitude, local_sidereal_time):
 n=-math.cos(local_sidereal_time)
 d=math.sin(obliquity_of_ecliptic)*math.tan(observer_latitude) + math.cos(obliquity_of_ecliptic)*math.sin(local_sidereal_time)
 return math.atan(n/d)

#Page 99
def angle_between_ecliptic_and_horizon(obliquity_of_ecliptic, observer_latitude, local_sidereal_time):
 n=math.cos(obliquity_of_ecliptic)*math.sin(observer_latitude) - math.sin(obliquity_of_ecliptic)*math.cos(observer_latitude)*math.sin(local_sidereal_time)
 return math.acos(n)

#Page 100
def angle_between_ecliptic_and_north_celestial_pole(obliquity_of_ecliptic, ecliptical_latitude_of_star, ecliptical_longitude_of_star):
 n=math.cos(ecliptical_longitude_of_star)*math.tan(obliquity_of_ecliptic)
 d=math.sin(ecliptical_latitude_of_star)*math.sin(ecliptical_longtitude_of_star)*math.tan(obliquity_of_ecliptic)-math.cos(ecliptical_latitude_of_star)
 return math.atan(n/d)

#by setting the ecliptical latitude to zero, we get the angle between the ecliptic and and the east-west direction on the celestial sphere
def angle_between_ecliptic_and_celestial_sphere(ecliptical_longitude, obliquity_of_ecliptic):
 n = -math.cos(ecliptical_longitude)*math.tan(obliquity_of_ecliptic)
 return math.atan(n)

#the diurnal path of a body relative to the horizon
#this function neglects refraction, and the change in declination of the object.
#it can be off by up to 4' for the moon, and 1 degree for the sun
def dirunal_path_horizon(declination_of_body, observer_latitude):
 b=math.tan(declination_of_body)*math.tan(observer_latitude)
 c=math.sqrt(1-b**2)
 n=c*math.cos(declination_of_body)
 d=math.tan(observer_latitude)
 j=math.atan(n/d)
 return j

############
#Chapter 15#
############

#Page 101 & 102
#The value of 34' is used as the typical value of refraction

#This function does not take into account variations in refraction, and refers to the geometric center of the object
#Bodies for which a disk can be resolved are usually measured from their upper limb.
#So for the Sun, the semidiameter of 16' is added to the geometric center
#Variations in refraction can alter the rise and set times by around 20 seconds.
#For objects in orbit around the Earth, which can be resolved (like the ISS) some well documented decisions will need to be
#made about how where the center and limbs of the object are.

def hour_angle_of_rise_or_set_low_fidelity(observer_latitude,object_declination):
 n=-math.tan(observer_latitude)*math.tan(object_declination)
 return math.acos(n)

#Page 102 & 103

#This is the check Meeus provides in "note 2 at the end of this chapter"
#Possible Bug Alert#
#Meeus says to check the second member of the function.  It is not clear if he means the second term in the numerator,
#or the second term of the equation if you rewrite it as the sum of two fractions.  I'll find out in testing
def is_circumpolar(observer_latitude,object_declination):
 if abs(math.sin(observer_latitude)*math.sin(object_declination)) > 1:
  return True
 else:
  return False

#sidereal_time is 0h Universal Time on the day D of interest
#ra_ and dec_ are the right ascension and declination of the object at 0h Dynamical Time for D-1, D, and D+1 respectively
def rise_transit_set(sidereal_time,ra1,ra2,ra3,dec1,dec2,dec3):
 #first we check to see if the object is circumpolar.  If it is, it will not have a rise or set time.
 #next we calculate the approximate times
 return 'nothing yet'

#################################################################################
#These are functions I wrote to make processing easier. They are not from Meeus.#
#################################################################################

def angle_convert(angle, c_from, c_to):
 if c_from.lower() == 'dms':
  if c_to.lower() == 'dd':
   return angle[0]+(angle[1]+(angle[2]/60.0))/60.0
  elif c_to.lower() == 'rad':
   return math.pi*angle_convert(angle,'dms','dd')/180.0
  else:
   return None
 if c_from.lower() == 'dd':
  if c_to.lower() == 'rad':
   return math.pi*angle/180.0
  elif c_to.lower() == 'dms':
   a = [1,2,3]	#cheap init of the list
   a[0]=int(angle)
   b=(angle-a[0])/60.0
   a[1]=int(b)
   c=(b-a[1])/60.0
   return [a,b,c]
  else:
   return None
 if c_from.ower() == 'rad':
  if c_to.lower() == 'dd':
   return 180.0*angle/math.pi
  elif c_to.lower() == 'dms':
   angle_convert(angle_convert(angle,'rad','dd'),'dd','dms')
  else:
   return None

def add_dms(angle1, angle2):
 c=[0,0,0]
 for i in range(3):
  c[i]=angle1[i]+angle2[i]
 if i>0 and c[i]>60:
  c[i]=c[i]-60
  c[i-1]+=1
 return c

#Returns True id the year is a leap year
def is_leap_year(year):
 if year%4 > 0:
  return False
 elif year%100 > 0:
  return True
 elif year%400 > 0:
  return True

#Converts a decimal day fraction to [hours,minutes,seconds]
def hours_minutes_seconds(day_fraction):
 hours = int(day_fraction*24)
 minutes = int((day_fraction*24 - hours)*60)
 seconds = ((day_fraction*24 - hours)*60 - minutes)*60
 return [hours, minutes, seconds]

#Returns the day fraction equivalent to the given time of day
def day_fraction(hours, minutes,seconds):
 a=minutes+seconds/60
 b=hours+a/60
 return b/24

#This funtion takes the month's name, abbreviation, or number and converts it to the other
#Output formats are 'name, 'number', and 'abbreviation
#If the function receives invalid input, it will raise an exception. This is not considered a bug. It will help debug the code that called the function
def month_name_number(month, output):
 names = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'November', 'December']
 abbreviations = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
 if output == 'name':
  try:
   return names[month-1]
  except:
   m=month_name_number(month,'number')
   return names[m-1]
 elif output == 'abbreviation':
  try:
   return abbreviations[month-1]
  except:
   m=month_name_number(month,'number')
 elif output == 'number':
  try:
   n=[i for (i,x) in enumerate(names) if x.lower() == month.lower()] or [i for (i,x) in enumerate(abbreviates) if x.lower() == month.lower()]
   return n[0]+1
  except:
   return month
  else:
   return None