
//set dirichlet nodes to all voxels outside of a radius, and remove any outside that radius
void activate_nodes(double radius_of_activation, std::vector<double> center_of_activation, std::vector<double> new_values)
{
  if(new_values.size()!=microenvironment.number_of_densities())
  {
    std::cout<<"WARNING!! Mismatched number of density values! \n";
    return;
  }
	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )//
	{
		if(norm(microenvironment.voxels(i).center, center_of_activation)>=radius_of_activation)
		{
		  microenvironment.add_dirichlet_node(i, new_values);//<<-- is this right?
      microenvironment.density_vector(i)=new_values;
		}
		else
		{
			microenvironment.remove_dirichlet_node(i);
		}
	}
  return;
}
