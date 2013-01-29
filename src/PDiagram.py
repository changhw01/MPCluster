import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


class PDiagram:

    def __init__(self,cM,sM,rep_clst):    

        ''' Initialization

            cM: cluster version
            sM: size of cluster 
            rep_clst: the representation of cluster 
            
        '''

        self.cM = cM
        self.sM = sM
        self.rep_clst = rep_clst

        self.n_h, self.n_r = cM.shape

    def show(self,ax=None):
        ''' Draw the persistence diagram 
        
            ax: the axis object 
        '''
    
        # prepare for plotting
        if ax is None:
            fig = plt.figure()
            cmap = plt.cm.summer
            ax = fig.add_subplot(111,aspect='equal')        

        # base coordinates
        centers_x = np.array(range(1,self.n_r+1))
        centers_y = np.array(range(self.n_h,-1,-1))

        gap_conn = 0.2
        gap_disc = 0.4
        base_rec_width = gap_disc*2
        base_rec_height = gap_disc*2

        base_coords_x_left = centers_x-gap_disc
        base_coords_x_right = centers_x+gap_disc
        base_coords_y_bottom = centers_y-gap_disc
        base_coords_y_up = centers_y+gap_disc


        # for color
        colorM = self.sM/ float(np.max(self.sM))

        # one rectangle for each clustering 
        for i in range(self.n_h):
            for j in range(self.n_r):
                if self.cM[i,j]>0:
                    x_left = base_coords_x_left[j]
                    y_bottom = base_coords_y_bottom[i]
            
                    rect = Rectangle((x_left,y_bottom),width=base_rec_width, height=base_rec_height,color=cmap(colorM[i,j]))
                    rect.set_linewidth(0)
                    ax.add_artist(rect)

        # fill up the vertical gaps if necessary
        for i in range(self.n_h):
            for j in range(1,self.n_r):
                if self.cM[i,j]==-1:
                    continue

                if self.cM[i,j-1]==self.cM[i,j]:
                    x_left = base_coords_x_right[j-1]
                    y_bottom = base_coords_y_bottom[i]

                    rect = Rectangle((x_left,y_bottom),width=gap_conn,height=base_rec_height,color=cmap(colorM[i,j]))
                    rect.set_linewidth(0)
                    ax.add_artist(rect)

        # fill up the horizonal gaps if necessary
        # pay attention for a 2x2 block
        for i in range(1,self.n_h):
            for j in range(self.n_r):

                if self.cM[i,j]==-1:
                    continue
                
                if self.cM[i-1,j]==self.cM[i,j]:
                    x_left = base_coords_x_left[j]
                    y_bottom = base_coords_y_up[i]

                    if (j<(self.n_r-1)) and (self.cM[i,j]==self.cM[i,j+1]) and (self.cM[i,j]==self.cM[i-1,j+1]):
                        rec_hlen = base_rec_width + gap_conn
                    else:
                        rec_hlen = base_rec_width

                    rect = Rectangle((x_left,y_bottom),width=rec_hlen,height=gap_conn,color=cmap(colorM[i,j]))
                    rect.set_linewidth(0)
                    ax.add_artist(rect)

        ax.autoscale_view()
        #ax.figure.canvas.draw()
        plt.axis([0, self.n_r+1, 0, self.n_h+1])
        fig.show()
        return ax
