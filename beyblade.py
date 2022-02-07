import RayTracer
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import glob
from PIL import Image
import os
import moviepy.editor as mp

def main():

    a=0.4
    initial_theta=-0.5
    for frame in range(180):
        initial_theta=initial_theta+0.5
        RayTracer.image(a,initial_theta,frame,wavelength_function=lambda_from_xy_3,unobservable_to_gray=True,set_intensity=False)
    
    gr_gif_frames=[]
    newtonian_gif_frames=[]

    gr_images=sorted(glob.glob("*gr_frame*"),key=os.path.getmtime)
    for i in gr_images:

        new_frame=Image.open(i)
       gr_gif_frames.append(new_frame)
    

    newtonian_images=sorted(glob.glob("*newtonian_frame*"),key=os.path.getmtime)
    for i in newtonian_images:

        new_frame=Image.open(i)
        newtonian_gif_frames.append(new_frame)
    

    gr_gif_frames[0].save('gr.gif',format='GIF',append_images=gr_gif_frames[1:],save_all=True,duration=100)
    newtonian_gif_frames[0].save('newtonian.gif',format='GIF',append_images=newtonian_gif_frames[1:],save_all=True,duration=100)
    clip=mp.VideoFileClip("gr.gif")
    clip.write_videofile("gr_visible.mp4")

    clip=mp.VideoFileClip("newtonian.gif")
    clip.write_videofile("newtonian_visible.mp4")


def lambda_from_xy_4(x,y):
    
    #beyblade
    a=0.4
    initial_lambda=550
    r=np.sqrt(x**2+y**2-a**2)    
    d_lambda=250
    spiral_a=0
    spiral_b=18/np.pi
    phi=np.arctan(y/x)
    if type(phi)!=type(np.zeros(1)):
        if phi<0:
            phi=phi+np.pi
    else:
        for i,placeholder in enumerate(phi):
            if placeholder<0:
                phi[i]=placeholder+np.pi
    spiral_par=r-(spiral_a+spiral_b*phi)
    true_lambda=initial_lambda-d_lambda*np.exp(-abs(spiral_par)/2.5)
    return true_lambda

if __name__=="__main__":
    main()