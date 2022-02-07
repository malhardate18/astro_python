import os
import platform
import sys
import numpy as np
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d as ax3d
from matplotlib import animation
import matplotlib.colors
#from moviepy.video.io.bindings import mplfig_to_npimage #this is because matplotlib acts up whenever you try to do a single thing


def show_progress_bar(cells,size,preText="",postText=""):

    #progress bar size
    n_bar=500
    cells=cells/size
    print([f" {preText} {'=' * int(n_bar * q):{n_bar}s}] {postText} "])
    return None

def compile_raytracer(cpp_file="raytracing.cpp",exe_file="raytracing.o"):
    
    #this function compiles the raytracing c++ file
    os.system("echo Compiling "+cpp_file)
    os.system("g++ " +cpp_file+" -o "+exe_file)
    print("Done.")
    return exe_file

def execute_raytracer(parameters,method="Simple",exe_file="raytracing.o",cpp_file="raytracing.cpp",not_compiled=False,time_get=False,j_get=False):
    #you can change the method in examples.py
    #in general,check if the C++ file is compiled. If not,compile it. Depending on your operating system,the way you run an executable changes
    #AIX OS not supported
    #the "parameters" input should be a 1x4 array
    #the first parameter is the x coordinate in the image plane
    #the second parameter is the y coordinate in the image plane
    #the third parameter is the inclination angle of the image plane
    #the fourth parameter is the normalized black hole spin a,satisfying |a|<1
   
    if not_compiled:
        exe_file=compile_raytracer(cpp_file,exe_file)
    elif os.path.isfile(exe_file):
        pass
    else:
        exe_file=compile_raytracer(cpp_file,exe_file)
    
    if platform.system()=="Windows":
        exe_file=".\\"+exe_file
    else:
        exe_file="./"+exe_file #linux or macOS
    
    if len(parameters)!=4:
        print("Got %d parameters (need 4)"%(len(parameters))) 
        raise ValueError("Got %d parameters (need 4)"%(len(parameters))) 
    
    
    parameters.append(method)

    allowed_methods=["Simple","Disk","NoDisk"]
    if parameters[-1] not in allowed_methods:
        print("Unallowed Method")
        print("Got %s method but only %s ,%s ,%s allowed."%(method,*allowed_methods)) 
        print("The method is now Simple (default)")
        parameters[-1]="Simple"
        
    parameters=[str(param) for param in parameters]
    
    pipe=os.popen(exe_file+" "+" ".join(parameters))#replace with subprocess?
    
    if time_get:
        if j_get:
            x,y,z,t,j=convert_cpp_output_to_numpy(pipe,a=float(parameters[-2]),time_get=time_get,j_get=j_get)
            return x,y,z,t,j
        x,y,z,t=convert_cpp_output_to_numpy(pipe,a=float(parameters[-2]),time_get=time_get,j_get=j_get)
        return x,y,z,t
    if j_get:
        x,y,z,j=convert_cpp_output_to_numpy(pipe,a=float(parameters[-2]),time_get=time_get,j_get=j_get)
        return x,y,z,j
    x,y,z=convert_cpp_output_to_numpy(pipe,a=float(parameters[-2]))
    return x,y,z
    
def convert_cpp_output_string_to_array(line):
    
    #hacky code that converts the output of the C++ backend into a numpy array

    arr=[]
    i=0
    for j in range(len(line)):
        if line[j]==",":
            arr.append(float(line[i:j]))
            i=j+1
    return np.array(arr)

def convert_cpp_output_to_numpy(pipe,a,time_get=False,j_get=False):
    
    #does the same thing but calls the function defined just above to convert things to numpy

    r_line=pipe.readline()
    theta_line=pipe.readline()
    phi_line=pipe.readline()
    t_line=pipe.readline()
    j_line=pipe.readline()
    
    r=convert_cpp_output_string_to_array(r_line)
    theta=convert_cpp_output_string_to_array(theta_line)
    phi=convert_cpp_output_string_to_array(phi_line) 
    
    if time_get:
        t=convert_cpp_output_string_to_array(t_line)
    
    if j_get:
        j=float(j_line)
    
    
    x,y,z=spherical_to_cartesian(r,theta,phi,a)
    
    if time_get:
        if j_get:
            return x,y,z,t,j
        return x,y,z,t
    if j_get:
        return x,y,z,j
    return x,y,z

#DO NOT TOUCH THESE TWO FUNCTIONS-PROGRAM MISBEHAVES IF EDITED
def spherical_to_cartesian(spherical_r,spherical_theta,spherical_phi,a):
    
    #self-explanatory
    cartesian_x=np.sqrt(a**2+spherical_r**2)*np.sin(spherical_theta)*np.cos(spherical_phi)
    cartesian_y=np.sqrt(a**2+spherical_r**2)*np.sin(spherical_theta)*np.sin(spherical_phi)
    cartesian_z=spherical_r*np.cos(spherical_theta)
    return cartesian_x,cartesian_y,cartesian_z

def cartesian_to_spherical(cartesian_x,cartesian_y,cartesian_z,a):
    
    #self-explanatory
    spherical_phi=np.arctan(cartesian_y/cartesian_x)
    convenient_R=np.sqrt(cartesian_x**2+cartesian_y**2+cartesian_z**2)
    spherical_r=np.sqrt(1/2*(convenient_R**2-a**2+np.sqrt((convenient_R**2-a**2)**2+4*a**2*cartesian_z**2)))
    spherical_theta=np.arccos(cartesian_z/spherical_r)
    return spherical_r,spherical_theta,spherical_phi

def compute_innermost_stable_circular_orbit(a):
	
    #gets the innermost stable circular orbit
	Z_1=1+(1-a**2)**(1/3)*((1+a)**(1/3)+(1-a)**(1/3)) 
	Z_2=np.sqrt(3*a**2+Z_1**2)
	return (3+Z_2-np.sign(a)*np.sqrt((3-Z_1)*(3+Z_1+2*Z_2)))

def compute_event_horizon(a):

    #gets the event horizon radius
    return 1+np.sqrt(1-a*a)

def prepare_to_draw(a,fig_width=5,fig_height=4,view_theta=30,view_phi=-90,axis="off"):
   
    #starts to prepare the final image - a black disk is drawn on a matplotlib figure
    fig=plt.figure(figsize=(fig_width,fig_height))
    ax=fig.add_subplot(111,projection="3d")
    ax.set_xlim(-15,15)
    ax.set_ylim(-15,15)
    ax.set_zlim(-10,10)
    ax.view_init(90-view_theta,view_phi)
    
    ax=draw_black_hole(a,ax)
    ax=draw_Disk(a,ax)
    
    plt.axis(axis)
    
    return fig,ax

def prepare_to_draw_no_Disk(a,fig_width=5,fig_height=4,view_theta=30,view_phi=-90,axis="off"):
    
    #same thing as above except the disk isn"t drawn
    fig=plt.figure(figsize=(fig_width,fig_height))
    ax=fig.add_subplot(111,projection="3d")
    ax.set_xlim(-15,15)
    ax.set_ylim(-15,15)
    ax.set_zlim(-9,9)
    ax.view_init(90-view_theta,view_phi)
    
    ax=draw_black_hole(a,ax)
    
    plt.axis(axis)
    
    return fig,ax


def draw_black_hole(a,ax):
    
    #just draw the black hole
    rH=compute_event_horizon(a)
    
    u,v=np.mgrid[0:2*np.pi:20j,0:np.pi:20j] #this is the linspace for the number of steps to iterate through except it isn"t actually a linspace
    xs=np.sqrt(rH*rH+a*a)*np.cos(u)*np.sin(v)
    ys=np.sqrt(rH*rH+a*a)*np.sin(u)*np.sin(v)
    zs=rH*np.cos(v)

    ax.plot_surface(xs,ys,zs,color="black",alpha=0.5) #alpha is the opacity of the black hole
    plt.axis("off")    
    
    return ax
    

def draw_Disk(a,ax):
   
    #this draws the disk
    #p is the innermost radius
    #q is the outermost radius
    innermost_stable_circular_orbit_radius=compute_innermost_stable_circular_orbit(a)
    
    N=100

    theta_d=np.linspace(0,2.*np.pi,N)
    phi_d=np.linspace(0,2.*np.pi,N)
    theta_d,phi_d=np.meshgrid(theta_d,phi_d)
    innermost_radius,outermost_radius=20,innermost_stable_circular_orbit_radius
    average_radius,mean_difference=(innermost_radius+outermost_radius)/2,(innermost_radius-outermost_radius)/2
    x_d=(average_radius+mean_difference*np.cos(theta_d))*np.cos(phi_d)
    y_d=(average_radius+mean_difference*np.cos(theta_d))*np.sin(phi_d)
    z_d=0.1*mean_difference*np.sin(theta_d)

    ax.contourf(x_d,y_d,z_d,[0.000000001,2],zdir="z",cmap=cm.inferno,alpha=0.5)
    
    return ax
    
def draw_photon_path(ax,x,y,z,color=None):
    
    #this function draws the photon path
    ax.plot(x,y,z,color="xkcd:sky blue")
    return


def start_animating_rays(xs,ys,zs,ts,a,cmap="jet",n_frame=None,Disk=False,burst_mode=False,lw=2,ls="-",interval=1,blit=False,repeat=True, fig_width=9,fig_height=6,view_theta=30,view_phi=-90):
    
    #this function animates ray paths, gettings values from x_i(t_i), y_i(t_i), z_i(t_i)
    #the input is lists, xs=[x1(t1)...] and so on
    
    if n_frame is None:
        n_frame=int(max([len(x) for x in xs])/3)
        
    total_number_of_rays=len(xs)
    
    if cmap is None:
        colors=[None for i in range(total_number_of_rays)]
    else:
        cm=plt.cm.get_cmap(cmap)
        colors=cm(np.linspace(0,1,total_number_of_rays))

    def get_t_i(x,y,z,t,initial_r=30):
        
        #this function gives you the starting time, i.e. when r=initial_r
        r,_,_=cartesian_to_spherical(x,y,z,a)
        ind=np.argmin(r>initial_r) # gets first instance of r<initial_r in array.
        t_i=t[ind]
        return t_i,ind
    
    def get_t_max(ts):
        
        #this function returns the largest proper time of all the paths
        tm=max([t[-1] for t in ts])
        return tm
        
    x_animate=[[] for i in range(total_number_of_rays)]
    y_animate=[[] for i in range(total_number_of_rays)]
    z_animate=[[] for i in range(total_number_of_rays)]
    
    lamb_animate=np.linspace(0,1,n_frame)
    
    t_j=get_t_max(ts)
    
    #this is where we interpolate
    for i in range(total_number_of_rays):
        t_i,ind=get_t_i(xs[i],ys[i],zs[i],ts[i])
        t_m_i=max(ts[i])
        if t_m_i!=t_j:
            interpolated_x=interp1d((np.append(ts[i][ind:],t_j)-t_i)/(t_j-t_i),np.append(xs[i][ind:],xs[i][-1]),kind="slinear") #duplicate error if anything higher is taken
            interpolated_y=interp1d((np.append(ts[i][ind:],t_j)-t_i)/(t_j-t_i),np.append(ys[i][ind:],ys[i][-1]),kind="slinear")
            interpolated_z=interp1d((np.append(ts[i][ind:],t_j)-t_i)/(t_j-t_i),np.append(zs[i][ind:],zs[i][-1]),kind="slinear")
        else:
            interpolated_x=interp1d((ts[i][ind:]-t_i)/(t_j-t_i),xs[i][ind:],kind="slinear")
            interpolated_y=interp1d((ts[i][ind:]-t_i)/(t_j-t_i),ys[i][ind:],kind="slinear")
            interpolated_z=interp1d((ts[i][ind:]-t_i)/(t_j-t_i),zs[i][ind:],kind="slinear")

        x_animate[i]=interpolated_x(lamb_animate)
        y_animate[i]=interpolated_y(lamb_animate)
        z_animate[i]=interpolated_z(lamb_animate)


    
    if Disk:
        fig,ax=prepare_to_draw(a,fig_width=fig_width,fig_height=fig_height,view_theta=view_theta,view_phi=view_phi)
    else:
        fig,ax=prepare_to_draw_no_Disk(a,fig_width=fig_width,fig_height=fig_height,view_theta=view_theta,view_phi=view_phi)
    
    line,=ax.plot([],[])
    plotlays,plotcols=[total_number_of_rays],colors

    lines=[]
    for index in range(total_number_of_rays):
        lobj=ax.plot([],[],lw=lw,color=plotcols[index],linestyle=ls)[0]
        lines.append(lobj)

    def init():
        for line in lines:
            line.set_data(np.asarray([]),np.asarray([]))
            line.set_3d_properties(np.asarray([]))
        return lines

    def animate(i):
        xlist=[[] for i in range(total_number_of_rays)]
        ylist=[[] for i in range(total_number_of_rays)]
        zlist=[[] for i in range(total_number_of_rays)]
        if burst_mode:
            for j in range(total_number_of_rays):
                xlist[j]=x_animate[j][max(0,i-int(n_frame/10)):i]
                ylist[j]=y_animate[j][max(0,i-int(n_frame/10)):i]
                zlist[j]=z_animate[j][max(0,i-int(n_frame/10)):i]
        else:
            for j in range(total_number_of_rays):
                xlist[j]=x_animate[j][:i]
                ylist[j]=y_animate[j][:i]
                zlist[j]=z_animate[j][:i]

        for lnum,line in enumerate(lines):
            line.set_data(np.asarray(xlist[lnum]),np.asarray(ylist[lnum]))
            line.set_3d_properties(np.asarray(zlist[lnum]))
        return lines

    anim=animation.FuncAnimation(fig,animate,init_func=init,frames=n_frame,interval=interval,blit=blit,repeat=repeat)
    anim.save("animation.gif")
    
    return fig,ax,anim


def make_wavelengths_visible(wavelength,gamma=0.8,alpha=0.1,unobservable_to_gray=True):
    
    #this function converts a given wavelength of light to an approximate RGB value, with non-visible wavelengths to be set to gray
        
    wavelength=float(wavelength)
    if wavelength>=380 and wavelength<=750:

        A=1
    else:
        if unobservable_to_gray:
            A=alpha
        else:
            A=1
    if wavelength<380:

        if unobservable_to_gray:
            wavelength=379
        else:
            wavelength=380
    
    if wavelength>750:

        if unobservable_to_gray:
            wavelength=751
        else:
            wavelength=750
    
    if wavelength>=380 and wavelength<=440:

        attenuation=0.3+0.7*(wavelength-380)/(440-380)
        R=((-(wavelength-440)/(440-380))*attenuation)**gamma
        G=0
        B=(1*attenuation)**gamma
    elif wavelength>=440 and wavelength<=490:

        R=0
        G=((wavelength-440)/(490-440))**gamma
        B=1
    elif wavelength>=490 and wavelength<=510:

        R=0
        G=1
        B=(-(wavelength-510)/(510-490))**gamma
    elif wavelength>=510 and wavelength<=580:

        R=((wavelength-510)/(580-510))**gamma
        G=1
        B=0
    elif wavelength>=580 and wavelength<=645:

        R=1
        G=(-(wavelength-645)/(645-580))**gamma
        B=0
    elif wavelength>=645 and wavelength<=750:

        attenuation=0.3+0.7*(750-wavelength)/(750-645)
        R=(1*attenuation)**gamma
        G=0
        B=0
    else:

        R=0
        G=0
        B=0
    return (R,G,B,A)

def compute_spectral_color_map(unobservable_to_gray=True):
    
    #This function returns a spectral color map to use further
    clim=(350,780)
    norm=plt.Normalize(*clim)
    wl=np.arange(clim[0],clim[1]+1,2)
    colorlist=list(zip(norm(wl),[make_wavelengths_visible(w,unobservable_to_gray=unobservable_to_gray) for w in wl]))
    spectralmap=matplotlib.colors.LinearSegmentedColormap.from_list("spectrum",colorlist)
    return spectralmap

def image(a,theta0,frame_counter,radial_photons=50,angular_photons=200,r_out=20,rest_wavelength=550,wavelength_function=None,unobservable_to_gray=False,set_intensity=False,print_progress=False):
    #radial_photons*angular_photons=total photons
    
    #This function is the actual camera image used
    save_string=["radial_photons=",str(radial_photons),"angular_photons=",str(angular_photons),"a=",str(a),"theta=",str(theta0),"r_out=",str(r_out)]
    save_string="".join(save_string)
        
    files_exist=os.path.isfile("im_x"+save_string+".npy")
    files_exist=files_exist and os.path.isfile("im_y"+save_string+".npy")
    files_exist=files_exist and os.path.isfile("phys_y"+save_string+".npy")
    files_exist=files_exist and os.path.isfile("phys_y"+save_string+".npy")
    files_exist=files_exist and os.path.isfile("red_shift"+save_string+".npy")
    
    if files_exist:
        print("The raytracer is loading some physical parameters")
        ri=compute_innermost_stable_circular_orbit(a)
        rH=compute_event_horizon(a)
        
        im_x=np.load("im_x"+save_string+".npy")
        im_y=np.load("im_y"+save_string+".npy")
        phys_x=np.load("phys_x"+save_string+".npy")
        phys_y=np.load("phys_y"+save_string+".npy")
        red_shifts=np.load("red_shift"+save_string+".npy")
        
    else:
        print("Computing %d photons - hang on!"%(radial_photons*angular_photons))
        def U0(r,a):
            return (1+a*np.power(r,-3/2))/np.sqrt(1-3/r+2*a*np.power(r,-3/2))
        def Uphi(r,a):
            return 1/(np.power(r,3/2)*np.sqrt(1-3/r+2*a*np.power(r,-3/2)))
        def g(r,a,j,ri,r_out=20):
            if r<ri:
                return 10
            elif r>r_out:
                return 0
            return 1/(U0(r,a))*1/(1+j*(Uphi(r,a)/U0(r,a))) 
    
        ri=compute_innermost_stable_circular_orbit(a)
        rH=compute_event_horizon(a)
    
        im_r=[ri+(r_out+5-ri)*i/(radial_photons-1) for i in range(radial_photons)]
        im_phi=[np.pi/2+2*np.pi*i/(angular_photons-1) for i in range(angular_photons)]

        im_x=[]
        im_y=[]
        red_shifts=[]
        phys_x=[]
        phys_y=[]
        
        f_min=0.1
        f_max=20
        total_number_of_rays_traced=0
        for i in range(radial_photons):
            for j in range(angular_photons):
                total_number_of_rays_traced += 1
                alpha=im_r[i]*np.cos(im_phi[j])
                beta=im_r[i]*np.sin(im_phi[j])*np.cos(theta0*np.pi/180)
                if beta> 0:
                    beta *= np.power(1/np.cos(theta0*np.pi/180),.7)
                im_x.append(alpha)
                im_y.append(beta)
                x,y,z,j=execute_raytracer([alpha,beta,theta0,a],time_get=False,j_get=True)
                rf,_,_=cartesian_to_spherical(x[-1],y[-1],z[-1],a)#gives the final radial coordinate
                if rf<rH:
                    f=10
                elif z[-1]>0.5:
                    f=10
                else:
                    f=g(rf,a,j,ri,r_out)
                red_shifts.append(f)
                phys_x.append(x[-1])
                phys_y.append(y[-1])
                if print_progress:
                    if total_number_of_rays_traced % 100==0:
                        show_progress_bar(Q=total_number_of_rays_traced,size=radial_photons*angular_photons,preText=" - Currently raytracing - ",postText="Total number of photons traced: "+str(total_number_of_rays_traced)+"/"+str(radial_photons*angular_photons))
        
        if print_progress:            
            show_progress_bar(Q=radial_photons*angular_photons,size=radial_photons*angular_photons,preText=" - Currently raytracing - ",postText="Total number of photons traced: "+str(radial_photons*angular_photons)+"/"+str(radial_photons*angular_photons))
            print("")
        
        red_shifts=[red_shifts[i] if not (red_shifts[i]<f_min or red_shifts[i]> f_max) else np.nan for i in range(radial_photons*angular_photons)]
    
        save_string=["radial_photons=",str(radial_photons),"angular_photons=",str(angular_photons),"a=",str(a),"theta=",str(theta0),"r_out=",str(r_out)]
        save_string="".join(save_string)
        np.save("im_x"+save_string+".npy",im_x)
        np.save("im_y"+save_string+".npy",im_y)
        np.save("phys_x"+save_string+".npy",phys_x)
        np.save("phys_y"+save_string+".npy",phys_y)
        np.save("red_shift"+save_string+".npy",red_shifts) 
    
    print("Done. ")
    
    #plotting starts here
    if wavelength_function is not None:
        wavelengths=[wavelength_function(phys_x[i],phys_y[i])/red_shifts[i] if not np.isnan(red_shifts[i]) else np.nan for i in range(radial_photons*angular_photons)]
    else:
        wavelengths=[rest_wavelength/red_shifts[i] if not np.isnan(red_shifts[i]) else np.nan for i in range(radial_photons*angular_photons)]
    
    imx_lim_u=max([im_x[i] if not np.isnan(wavelengths[i]) else 0 for i in range(radial_photons*angular_photons)])
    imx_lim_l=min([im_x[i] if not np.isnan(wavelengths[i]) else 0 for i in range(radial_photons*angular_photons)])
    imy_lim_u=max([im_y[i] if not np.isnan(wavelengths[i]) else 0 for i in range(radial_photons*angular_photons)])
    imy_lim_l=min([im_y[i] if not np.isnan(wavelengths[i]) else 0 for i in range(radial_photons*angular_photons)])
        
    spectralmap=compute_spectral_color_map(unobservable_to_gray=unobservable_to_gray)
    
    fig=plt.figure(figsize=(16,9))
    ax=fig.add_subplot(111)
    
    if set_intensity:
        colors=spectralmap((np.array([wavelengths[i] if not np.isnan(wavelengths[i]) else 0 for i in range(radial_photons*angular_photons)])-350)/(780-350))
        colors[:,3]=[red_shifts[i]**3 if not np.isnan(red_shifts[i]) else 0.1 for i in range(radial_photons*angular_photons)]
        colors[:,3]=[.999*(colors[i,3]-min(colors[:,3])) /(max(colors[:,3])-min(colors[:,3]))  for i in range(radial_photons*angular_photons)]
        sc=ax.scatter(im_x,im_y,s=100,color=colors,edgecolors=None)
    else:
        sc=ax.scatter(im_x,im_y,s=100,c=wavelengths,cmap=spectralmap,vmin=349,vmax=781,edgecolors=None)
            
    ax.set_xlabel(r"$X$",fontsize=20)
    ax.set_ylabel(r"$Y$",fontsize=20,rotation=0)

    ax.set_xlim(imx_lim_l-6.5,imx_lim_u+6.5)
    ax.set_ylim(imy_lim_l-5,imy_lim_u+5)

    fig1=plt.figure(figsize=(16,9))
    ax1=fig1.add_subplot(111)

    ax1.set_xlabel(r"$X$",fontsize=20)
    ax1.set_ylabel(r"$Y$",fontsize=20,rotation=0)
    
    ax1.set_xlim(imx_lim_l-30,imx_lim_u+30)
    ax1.set_ylim(imy_lim_l-30,imy_lim_u+30)
    
    if wavelength_function is None:
        plot_phi=np.linspace(0,2*np.pi,100)

        Disk_x_in=ri*np.cos(plot_phi)
        Disk_y_in=ri*np.cos(theta0*np.pi/180)*np.sin(plot_phi)

        Disk_x_out=r_out*np.cos(plot_phi)
        Disk_y_out=r_out*np.cos(theta0*np.pi/180)*np.sin(plot_phi)

        ax1.plot(Disk_x_in,Disk_y_in,color=make_wavelengths_visible(rest_wavelength),lw=0)
        ax1.plot(Disk_x_out,Disk_y_out,color=make_wavelengths_visible(rest_wavelength),lw=0)
        ax1.fill(np.append(Disk_x_in,Disk_x_out[::-1]),np.append(Disk_y_in,Disk_y_out[::-1]),color=make_wavelengths_visible(rest_wavelength))

        circle1=plt.Circle((0,0),rH,color="k")

        ax1.add_artist(circle1)

        if rH>ri*np.cos(theta0*np.pi/180):
            n=50
            for i in range(n):
                plot_phi=np.linspace(np.pi/2,3*np.pi/2,100)
                r=ri+(rH/np.cos(theta0*np.pi/180)-ri)*k/(n-1)
                Disk_x_line=r*np.sin(plot_phi)
                Disk_y_line=r*np.cos(theta0*np.pi/180)*np.cos(plot_phi)
                ax1.plot(Disk_x_line,Disk_y_line,color=make_wavelengths_visible(rest_wavelength),lw=1)
    else:
        n=50
        for i in range(n):
            plot_phi=np.linspace(0,2*np.pi,200)
            r=ri+(r_out-ri)*i/(n-1)
            Disk_x_line=r*np.cos(plot_phi)
            Disk_y_line=r*np.cos(theta0*np.pi/180)*np.sin(plot_phi)
            ax1.scatter(Disk_x_line,Disk_y_line,s=100,c=wavelength_function(Disk_x_line,Disk_y_line/np.cos(theta0*np.pi/180)),cmap=spectralmap,vmin=350,vmax=780)
    
        if rH> ri*np.cos(theta0*np.pi/180):
            for i in range(n):
                plot_phi=np.linspace(0,np.pi,200)
                r=rH*i/(n-1)
                Disk_x_line=r*np.cos(plot_phi)
                Disk_y_line=r*np.sin(plot_phi)
                ax1.scatter(Disk_x_line,Disk_y_line,s=1,color="k")
        else:
            for i in range(n):
                plot_phi=np.linspace(0,2*np.pi,400)
                r=rH*i/(n-1)
                Disk_x_line=r*np.cos(plot_phi)
                Disk_y_line=r*np.sin(plot_phi)
                ax1.scatter(Disk_x_line,Disk_y_line,s=1,color="k")

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax1.xaxis.set_visible(False)
    ax1.yaxis.set_visible(False)
    for annoying_line in ["top","right","left","bottom"]:
        ax.spines[annoying_line].set_visible(False)

    gr_frame=f"gr_frame{frame_counter}.png"
    newtonian_frame=f"newtonian_frame{frame_counter}.png"
    fig.savefig(gr_frame,bbox_inches="tight",pad_inches=0)
    fig1.savefig(newtonian_frame,bbox_inches="tight",pad_inches=0)
    plt.close(fig)
    plt.close(fig1)


    ax.xaxis.set_visible(True)
    ax.yaxis.set_visible(True)
    ax1.xaxis.set_visible(True)
    ax1.yaxis.set_visible(True)

    for annoying_line in ["top","right","left","bottom"]:
        ax.spines[annoying_line].set_visible(True)

    ax.set_title("General Relativity",fontsize=15,family="monospace")
    ax1.set_title("Classical Mechanics",fontsize=15,family="monospace")
    
    return fig,ax,fig1,ax1


def draw_photon_paths_from_parameters(spins,thetas,alphas,betas,method="NoDisk",fig_width=9,fig_height=6,view_theta=30,view_phi=-90):
    
    #convert parameters to lists
    if type(spins)!=type([]):
        spins=[spins]
    if type(thetas)!=type([]):
        thetas=[thetas]
    if type(alphas)!=type([]):
        alphas=[alphas]
    if type(betas)!=type([]):
        betas=[betas]
    
    #check lengths
    max_len=max([len(spins),len(thetas),len(alphas),len(betas)])
    len_test=(len(spins)==max_len or len(spins)==1)
    len_test=len_test and (len(thetas)==max_len or len(thetas)==1)
    len_test=len_test and (len(alphas)==max_len or len(alphas)==1)
    len_test=len_test and (len(betas)==max_len or len(betas)==1)
    
    if not len_test:
        raise ValueError("Initial ray parameters must be lists of length 1 or equal length. You have provided lists with lengths %d,%d,%d,%d."
        %(len(spins),len(thetas),len(alphas),len(betas)))
    
    if len(spins)==1:
        spins=[spins[0] for i in range(max_len)]
    if len(thetas)==1:
        thetas=[thetas[0] for i in range(max_len)]
    if len(alphas)==1:
        alphas=[alphas[0] for i in range(max_len)]
    if len(betas)==1:
        betas=[betas[0] for i in range(max_len)]
    
    radial_photonsay=max_len
    xs,ys,zs,ts=[[] for i in range(radial_photonsay)],[[] for i in range(radial_photonsay)],[[] for i in range(radial_photonsay)],[[] for i in range(radial_photonsay)]

    for i,a0,b0,theta0,a in zip(range(radial_photonsay),alphas, betas,thetas,spins):
        x,y,z,t=execute_raytracer([a0,b0,theta0,a],method,time_get=True)
        xs[i],ys[i],zs[i],ts[i]=x,y,z,t
    
    if method=="NoDisk":
        fig,ax=prepare_to_draw_no_Disk(spins[0],fig_width=fig_width,fig_height=fig_height,view_theta=view_theta,view_phi=view_phi)
    else:
        fig,ax=prepare_to_draw(spins[0],fig_width=fig_width,fig_height=fig_height,view_theta=view_theta,view_phi=view_phi)
    
    cm=plt.cm.get_cmap("jet")
    colors=cm(np.linspace(0,1,radial_photonsay))

    for i,x,y,z in zip(range(radial_photonsay),xs,ys,zs):
        draw_photon_path(ax,x,y,z,color=colors[i])

    return fig,ax

def start_animating_rays_from_parameters(spins,thetas,alphas,betas,method="NoDisk",cmap="jet",n_frame=None,burst_mode=False,
                lw=2,ls="-",interval=1,blit=False,repeat=True,fig_width=9,fig_height=6,view_theta=30,view_phi=-90):
    
    #turn parameters into lists
    if type(spins)!=type([]):
        spins=[spins]
    if type(thetas)!=type([]):
        thetas=[thetas]
    if type(alphas)!=type([]):
        alphas=[alphas]
    if type(betas)!=type([]):
        betas=[betas]
    
    #check parameter lengths
    max_len=max([len(spins),len(thetas),len(alphas),len(betas)])
    len_test=(len(spins)==max_len or len(spins)==1)
    len_test=len_test and (len(thetas)==max_len or len(thetas)==1)
    len_test=len_test and (len(alphas)==max_len or len(alphas)==1)
    len_test=len_test and (len(betas)==max_len or len(betas)==1)
    
    if not len_test:
        raise ValueError("Your list lengths are incorrect (must be 1 or l). Current list lengths: %d,%d,%d,%d."%(len(spins),len(thetas),len(alphas),len(betas)))
    
    if len(spins)==1:
        spins=[spins[0] for i in range(max_len)]
    if len(thetas)==1:
        thetas=[thetas[0] for i in range(max_len)]
    if len(alphas)==1:
        alphas=[alphas[0] for i in range(max_len)]
    if len(betas)==1:
        betas=[betas[0] for i in range(max_len)]
    
    radial_photonsay=max_len
    xs,ys,zs,ts=[[] for i in range(radial_photonsay)],[[] for i in range(radial_photonsay)],[[] for i in range(radial_photonsay)],[[] for i in range(radial_photonsay)]

    for i,a0,b0,theta0,a in zip(range(radial_photonsay),alphas, betas,thetas,spins):
        x,y,z,t=execute_raytracer([a0,b0,theta0,a],method,time_get=True)
        xs[i],ys[i],zs[i],ts[i]=x,y,z,t
    
    if method=="NoDisk":
        anim_fig,anim_ax,anim=start_animating_rays(xs,ys,zs,ts,spins[0],Disk=False,cmap=cmap,n_frame=n_frame,burst_mode=burst_mode,lw=lw,ls=ls,interval=interval,blit=blit,repeat=repeat,fig_width=fig_width,fig_height=fig_height,view_theta=view_theta,view_phi=view_phi)
    else:
        anim_fig,anim_ax,anim=start_animating_rays(xs,ys,zs,ts,spins[0],Disk=True,cmap=cmap,n_frame=n_frame,burst_mode=burst_mode,lw=lw,ls=ls,interval=interval,blit=blit,repeat=repeat,fig_width=fig_width,fig_height=fig_height,view_theta=view_theta,view_phi=view_phi)

    return anim_fig,anim_ax,anim

def plot_and_start_animating_rays_from_parameters(spins,thetas,alphas,betas,method="NoDisk",cmap="jet",n_frame=None,burst_mode=False,lw=2,ls="-",interval=1,blit=False,repeat=True,fig_width=9,fig_height=6,view_theta=60,view_phi=-90):
    
    #make parameters into lists
    if type(spins)!=type([]):
        spins=[spins]
    if type(thetas)!=type([]):
        thetas=[thetas]
    if type(alphas)!=type([]):
        alphas=[alphas]
    if type(betas)!=type([]):
        betas=[betas]
    
    ##make sure they have the right lengths
    max_len=max([len(spins),len(thetas),len(alphas),len(betas)])
    len_test=(len(spins)==max_len or len(spins)==1)
    len_test=len_test and (len(thetas)==max_len or len(thetas)==1)
    len_test=len_test and (len(alphas)==max_len or len(alphas)==1)
    len_test=len_test and (len(betas)==max_len or len(betas)==1)
    
    if not len_test:
        raise ValueError("Your list lengths are incorrect (must be 1 or l). Current list lengths:%d,%d,%d,%d."%(len(spins),len(thetas),len(alphas),len(betas)))
    
    if len(spins)==1:
        spins=[spins[0] for i in range(max_len)]
    if len(thetas)==1:
        thetas=[thetas[0] for i in range(max_len)]
    if len(alphas)==1:
        alphas=[alphas[0] for i in range(max_len)]
    if len(betas)==1:
        betas=[betas[0] for i in range(max_len)]
    
    radial_photonsay=max_len
    xs,ys,zs,ts=[[] for i in range(radial_photonsay)],[[] for i in range(radial_photonsay)],[[] for i in range(radial_photonsay)],[[] for i in range(radial_photonsay)]

    for i,a0,b0,theta0,a in zip(range(radial_photonsay),alphas, betas,thetas,spins):
        x,y,z,t=execute_raytracer([a0,b0,theta0,a],method,time_get=True)
        xs[i],ys[i],zs[i],ts[i]=x,y,z,t
    
    if method=="NoDisk":
        plot_fig,plot_ax=prepare_to_draw_no_Disk(spins[0],fig_width=fig_width,fig_height=fig_height,view_theta=view_theta,view_phi=view_phi)
    else:
        plot_fig,plot_ax=prepare_to_draw(spins[0],fig_width=fig_width,fig_height=fig_height,view_theta=view_theta,view_phi=view_phi)
    
    cm=plt.cm.get_cmap("jet")
    colors=cm(np.linspace(0,1,radial_photonsay))

    for i,x,y,z in zip(range(radial_photonsay),xs,ys,zs):
        draw_photon_path(plot_ax,x,y,z,color=colors[i])
    
    if method=="NoDisk":
        anim_fig,anim_ax,anim=start_animating_rays(xs,ys,zs,ts,spins[0],Disk=False,cmap=cmap,n_frame=n_frame,burst_mode=burst_mode,lw=lw,ls=ls,interval=interval,blit=blit,repeat=repeat,fig_width=fig_width,fig_height=fig_height,view_theta=view_theta,view_phi=view_phi)
    else:
        anim_fig,anim_ax,anim=start_animating_rays(xs,ys,zs,ts,spins[0],Disk=True,cmap=cmap,n_frame=n_frame,burst_mode=burst_mode,lw=lw,ls=ls,interval=interval,blit=blit,repeat=repeat,fig_width=fig_width,fig_height=fig_height,view_theta=view_theta,view_phi=view_phi)
        
    return plot_fig,plot_ax,anim_fig,anim_ax,anim



