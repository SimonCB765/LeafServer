import adjlistcreation
import Leafcull

def main(similarities, thresholdPercentage):
    """Perform the PDB redundancy removal.

    @param similarities: A record of the percentage sequence identity between the chains/entries up for culling.
    @type similarities : dictionary
    @param thresholdPercentage: The maximum permissible percentage sequence identity that any two chains/entries may possess.
    @type thresholdPercentage : float
    @param verboseOutput: Whether status updates should be printed out to the user.
    @type verboseOutput:  boolean
    
    """

    # Create the sparsematrix of the protein similarity graph.
    adjacent, proteinNames = adjlistcreation.main(similarities, thresholdPercentage)
    
    # Choose which proteins to remove from the similarity graph.
    if proteinNames == []:
        # This is True if there are no similarities greater than the given percentage sequence identity. If there are no
        # chains that are too similar, then there is no need to cull any chains from the network.
        proteinsToCull = []
    else:
        # Choose which chains to remove from the similarity graph.
        proteinsToCull, proteinsToKeep = Leafcull.main(adjacent, proteinNames)

    return proteinsToCull